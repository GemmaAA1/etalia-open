# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db import models
from django.conf import settings
from django.contrib.auth.models import BaseUserManager
from django.db.models import Q
from django.contrib.postgres.fields import ArrayField
from django.db.models.expressions import RawSQL

from gensim import matutils

from paperstream.core.models import TimeStampedModel
from paperstream.core.utils import pad_vector
from paperstream.altmetric.models import AltmetricModel
from paperstream.nlp.models import Model
from config.celery import celery_app as app
from .constants import FEED_STATUS_CHOICES, STREAM_METHODS_MAP
from .utils import *


logger = logging.getLogger(__name__)


class StreamManager(BaseUserManager):

    def create(self, **kwargs):
        """Create new UserFeed / UserFeedVector
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')

        papers_seed = None
        if 'papers_seed' in kwargs:
            papers_seed = kwargs.pop('papers_seed')

        obj = super(StreamManager, self).create(**kwargs)

        if papers_seed:
            obj.add_papers_seed(papers_seed)
            obj.save()

        return obj

    def create_main(self, **kwargs):
        """Create a new default UserFeed / UserFeedVector, named 'main' """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')
        return self.create(name='main', **kwargs)


class Stream(TimeStampedModel):
    """User Stream model

    NB: Design to have multiple feed per user,
    """

    # stream name
    name = models.CharField(max_length=100, default='main')

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='streams')

    # Papers that are used in the similarity matching
    seeds = models.ManyToManyField(Paper, through='StreamSeeds',
                                   related_name='stream_seeds')

    # relevant matches matched
    matches = models.ManyToManyField(Paper, through='StreamMatches',
                                     related_name='stream_matches')

    state = models.CharField(max_length=3, blank=True, default='NON',
                             choices=FEED_STATUS_CHOICES)

    message = models.CharField(max_length=127, blank=True, default='')

    objects = StreamManager()

    class Meta:
        unique_together = (('name', 'user'),)

    def set_state(self, state):
        self.state = state
        if state == 'IDL':
            self.set_message('')
        self.save()

    def set_message(self, message):
        self.message = message
        self.save()

    def __str__(self):
        return '{stream}@{email}'.format(stream=self.name,
                                         email=self.user.email)

    def add_papers_seed(self, papers):
        """Add paper to seed in bulk"""
        if papers:
            objs = []
            seed_pks = self.seeds.all().values_list('pk', flat=True)
            for paper in papers:
                if paper.pk not in seed_pks:
                    objs.append(StreamSeeds(stream=self, paper=paper))
            StreamSeeds.objects.bulk_create(objs)

        # Update stream vector
        self.update_stream_vector()

    def update_stream_vector(self):
        """Update stream vector for each NLP model"""
        model_pks = Model.objects\
            .filter(is_active=True)\
            .values_list('pk', flat='True')
        for model_pk in model_pks:
            ufv, _ = self.vectors.get_or_create(model_id=model_pk)
            ufv.update_vector()

    def clear_old_papers(self):
        """Remove outdating matches based on user settings"""

        # get time lapse corresponding to user settings
        from_date = (timezone.now() -
                     timezone.timedelta(
                         days=self.user.settings.stream_time_lapse)).date()

        # clean old matches
        StreamMatches.objects\
            .annotate(date=RawSQL("SELECT LEAST(date_ep, date_fs, date_pp) "
                                      "FROM library_paper "
                                      "WHERE id = paper_id", []))\
            .filter(Q(stream=self) &
                    Q(date__lt=from_date))\
            .delete()

    def clear_all(self):
        """Delete all matched matches"""
        StreamMatches.objects.filter(stream=self).delete()

    def log(self, level, message, options):
        """Local wrapper around logger"""
        if level == 'info':
            logger.info('Stream ({pk}/{stream_name}@{user_email}): '
                        '{message} - {options}'.format(
                pk=self.id,
                stream_name=self.name,
                user_email=self.user.email,
                message=message,
                options=options))

        elif level == 'debug':
            logger.debug('Feed ({pk}/{stream_name}@{user_email}): '
                         '{message} - {options}'.format(
                pk=self.id,
                stream_name=self.name,
                user_email=self.user.email,
                message=message,
                options=options))
        else:
            raise ValueError('level unknown')

    def reset(self):
        """reset stream"""
        self.clear_all()
        self.update()

    def update(self, restrict_journal=False):
        """Update Stream

        - Get all matches in time range excluding untrusted and paper already
        in user lib
        - Score matches
        - Get only the N top scored matches
        - Create/Update StreamPaper
        """

        self.set_state('ING')
        self.log('info', 'Updating', 'starting...')
        self.set_message('Cleaning')

        self.clear_old_papers()

        # get paper that remain and need to be updated
        self.set_message('Fetching data')
        ufp_pks_to_update = self.matches.values_list('pk', flat=True)
        # get user library
        lib_pks = self.user.lib.papers.values_list('pk', flat=True)

        # Retrieve scoring class
        Score = eval(dict(STREAM_METHODS_MAP)[self.user.settings.stream_method])
        # and instantiate
        scoring = Score(
            model=self.user.settings.stream_model,
            user=self.user)

        # Retrieve matches to score
        self.log('debug', 'Updating', 'fetching relevant matches...')
        # get mostsimilar task
        try:
            ms_task = app.tasks['paperstream.nlp.tasks.mostsimilar_{name}'.format(
                name=self.user.settings.stream_model.name)]
        except KeyError:
            raise KeyError
        # get user journals
        if restrict_journal:
            journal_pks = self.user.lib.journals.all().values_list('pk', flat=True)
        else:
            journal_pks = None
        # get neighbors of seed matches
        seed_pks = self.seeds.all().values_list('pk', flat=True)
        res = ms_task.delay('get_partition',
                             paper_pks=seed_pks,
                             time_lapse=self.user.settings.stream_time_lapse,
                             k=settings.FEED_K_NEIGHBORS,
                             journal_pks=journal_pks)
        # wait for Results
        target_seed_pks = res.get()
        # build target matches list
        target_pks = list(ufp_pks_to_update) + \
                     [pk for pk in target_seed_pks
                      if pk not in list(lib_pks) + list(ufp_pks_to_update)]

        # compute scores
        objs_list = []
        if target_pks:
            # compute scores
            self.log('debug', 'Updating', 'preparing...')
            scoring.prepare(seed_pks, target_pks)
            self.log('debug', 'Updating', 'scoring...')
            pks, scores = scoring.score()
            # sort scores
            self.log('debug', 'Updating', 'sorting...')
            nb_papers = int(settings.FEED_SIZE_PER_DAY *
                            self.user.settings.stream_time_lapse)
            best = matutils.argsort(scores,
                                    topn=nb_papers,
                                    reverse=True)
            # reshape
            results = [(pks[ind], float(scores[ind])) for ind in best]
            self.log('debug', 'Updating', 'done')

            # create/update UserFeedPaper
            self.log('debug', 'Updating', 'populating...')
            for pk, val in results:
                # update
                if pk in ufp_pks_to_update:
                    ufp, _ = StreamMatches.objects.get_or_create(
                        stream=self,
                        paper_id=pk)
                    ufp.score = val
                    ufp.is_score_computed = True
                    ufp.save()
                else:  # create in bulk
                    objs_list.append(StreamMatches(
                        stream=self,
                        paper_id=pk,
                        score=val,
                        is_score_computed=True))
            # bulk create
            StreamMatches.objects.bulk_create(objs_list)
            self.log('debug', 'Updating', 'done')
        self.set_state('IDL')

        self.log('info', 'Updating', 'DONE')


class StreamSeeds(TimeStampedModel):
    """Relationship table between stream and seed matches"""

    stream = models.ForeignKey(Stream)

    paper = models.ForeignKey(Paper)

    class Meta:
        unique_together = ('stream', 'paper')

    def clean(self):
        """check if paper is in user.lib"""
        assert self.paper.pk in \
            self.stream.user.lib.papers.values_list('pk', flat=True)


class StreamMatches(TimeStampedModel):
    """Relationship table between stream and matched matches with scoring"""

    stream = models.ForeignKey(Stream)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    is_score_computed = models.BooleanField(default=False)

    class Meta:
        ordering = ['-score']
        unique_together = [('stream', 'paper')]

    def print_score(self):
        return '{0:.1f}'.format(self.score * 100)

    def __str__(self):
        return '{paper}/{score}'.format(paper=self.paper.short_title,
                                        score=self.score)


class StreamVector(TimeStampedModel):
    """Feature vector for stream is defined as the averaged of paper vectors in
    stream
    """
    stream = models.ForeignKey(Stream, related_name='vectors')

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True)

    class Meta:
        unique_together = ('stream', 'model')

    def __str__(self):
        return '{stream_name}/{model_name}'.format(stream_name=self.stream.name,
                                                   model_name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        if self.vector:
            return self.vector[:self.model.size]
        else:
            return None

    def update_vector(self):
        vectors = self.stream.seeds\
            .filter(vectors__model=self.model)\
            .values_list('vectors__vector', flat=True)

        # bag of vectors
        vector = [sum(x) for x in zip(*vectors) if x[0]]  # not considering None
        self.set_vector(vector)
        return vector[:self.model.size]


class Trend(TimeStampedModel):
    """Trend feed retrieving trending papers in a broaden neighborhood"""

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='trends')

    name = models.CharField(max_length=100, default='main')

    matches = models.ManyToManyField(Paper, through='TrendMatches')

    # top n closest matches
    top_n_closest = models.IntegerField(default=5000)

    # Add very high trending and recent papers
    # look for trendy altmetric in the past n days
    time_lapse_top_altmetric = models.IntegerField(default=15)
    score_threshold = models.FloatField(default=1000.0)

    def __str__(self):
        return self.user.email

    class Meta:
        unique_together = (('user', 'name'),)

    def clear(self):
        """Remove all paper relations"""
        TrendMatches.objects.filter(trend=self).all().delete()

    def update(self):

        # clear stream
        self.clear()

        # Get mostsimilar task
        try:
            ms_task = app.tasks['paperstream.nlp.tasks.mostsimilar_{name}'.format(
                name=self.user.settings.trend_model.name)]
        except KeyError:
            raise KeyError

        # get main stream signature vector
        stream = Stream.objects.get(user=self.user, name='main')
        sv, _ = StreamVector.objects\
            .get_or_create(stream=stream, model=self.user.settings.trend_model)
        vec = sv.update_vector()

        # Submit Celery Task to get k_neighbors from stream signature
        res = ms_task.delay('knn_search',
                            seed=vec,
                            time_lapse=self.user.settings.trend_time_lapse,
                            top_n=self.top_n_closest)
        target_seed_pks = res.get()

        # Get altmetric scores for target_seed_pks and top matches
        d = timezone.now().date() - \
            timezone.timedelta(days=self.time_lapse_top_altmetric)
        # get paper in the neighborhood of vec of lib ranked by alt score +
        # add very high impact recent altmetric matches
        nb_papers = int(settings.FEED_SIZE_PER_DAY * \
                        self.user.settings.trend_time_lapse)
        data_altm = AltmetricModel.objects\
            .annotate(date=RawSQL("SELECT LEAST(date_ep, date_fs, date_pp) "
                                  "FROM library_paper "
                                  "WHERE id = paper_id", []))\
            .filter(Q(paper_id__in=target_seed_pks) |
                    (Q(date__gt=d) & Q(score__gt=self.score_threshold)))\
            .order_by('-score')\
            .values('paper__pk', 'score')[:nb_papers]

        # Populate trendfeedpaper
        obj_list = []
        for p in data_altm:
            obj_list.append(TrendMatches(
                        trend=self,
                        paper_id=p['paper__pk'],
                        score=p['score']))
        TrendMatches.objects.bulk_create(obj_list)


class TrendMatches(TimeStampedModel):

    trend = models.ForeignKey(Trend)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.0)

    class Meta:
        unique_together = [('trend', 'paper'), ]

