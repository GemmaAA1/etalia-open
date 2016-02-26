# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db import models

from django.contrib.auth.models import BaseUserManager
from django.db.models import Q
from django.contrib.postgres.fields import ArrayField
from django.db.transaction import atomic

from paperstream.core.models import TimeStampedModel
from paperstream.core.utils import pad_vector
from paperstream.nlp.models import MostSimilar
from paperstream.last_seen.models import LastSeen
from .constants import FEED_STATUS_CHOICES, STREAM_METHODS_MAP, TREND_METHODS_MAP
from .scoring import *


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

    # last update of the stream
    last_update = models.DateTimeField(default=None, blank=True, null=True)

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

    def get_old_streammatches(self):
        """Remove outdating matches based on user settings"""

        # get time lapse corresponding to user settings
        from_date = (timezone.now() -
                     timezone.timedelta(
                         days=self.user.settings.stream_time_lapse)).date()

        # clean old matches
        return StreamMatches.objects.filter(stream=self, date__lt=from_date)

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

    def update(self, restrict_journal=False):
        """Update Stream

        - Get all matches in time range excluding untrusted and paper already
        in user lib
        - Score
        - Get only the N top scored matches
        - Create/Update StreamPaper
        """

        self.set_state('ING')
        self.log('info', 'Updating', 'starting...')
        self.set_message('Cleaning')

        # get out-dated StreamMatches
        del_objs = self.get_old_streammatches()

        # get paper_id list for StreamMatches that remain
        current_paper_list = self.streammatches_set\
            .exclude(id__in=del_objs)\
            .values_list('paper_id', flat=True)

        # Instantiate Score
        Score = eval(dict(STREAM_METHODS_MAP)[self.user.settings.stream_method])
        # and instantiate
        journal_ratio = MostSimilar.objects.get(is_active=True).journal_ratio
        # method_arg = self.user.settings.stream_method_args or {}
        method_arg = {
            'vector_weight': self.user.settings.stream_vector_weight,
            'author_weight': self.user.settings.stream_author_weight,
            'journal_weight': self.user.settings.stream_journal_weight,
        }
        scoring = Score(stream=self, journal_ratio=journal_ratio, **method_arg)

        # Score
        results, date = scoring.score()

        # create/update UserFeedPaper
        create_objs = []
        save_objs = []
        for pk, val in results:
            if pk in current_paper_list:
                sm = StreamMatches.objects.get(
                    stream=self,
                    paper_id=pk)
                sm.score = val
                save_objs.append(sm)
            else:
                create_objs.append(StreamMatches(
                    stream=self,
                    paper_id=pk,
                    score=val,
                    date=date[pk],
                    is_score_computed=True))

        # Get last user visit
        try:
            last_seen = LastSeen.objects.when(user=self.user)
        except LastSeen.DoesNotExist:
            last_seen = None
            pass

        # atomic update
        self.atomic_update(del_objs, save_objs, create_objs, last_seen=last_seen)

        self.set_state('IDL')
        self.log('info', 'Updating', 'DONE')

    @atomic
    def atomic_update(self, delete_objs, save_objs, create_objs, last_seen=None):
        if delete_objs:
            delete_objs.delete()
        if save_objs:
            for obj in save_objs:
                obj.save()
        if create_objs:
            StreamMatches.objects.bulk_create(create_objs)
        if last_seen:
            StreamMatches.objects.filter(created__lt=last_seen)\
                .update(new=False)


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

    date = models.DateField()

    is_score_computed = models.BooleanField(default=False)

    new = models.BooleanField(default=True)

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

    state = models.CharField(max_length=3, blank=True, choices=FEED_STATUS_CHOICES)

    def __str__(self):
        return self.user.email

    class Meta:
        unique_together = (('user', 'name'),)

    def set_state(self, state):
        self.state = state
        self.save()

    def get_old_trendmatches(self):
        """Remove outdating matches based on user settings"""

        # get time lapse corresponding to user settings
        from_date = (timezone.now() -
                     timezone.timedelta(
                         days=self.user.settings.stream_time_lapse)).date()

        # clean old matches
        return TrendMatches.objects\
            .filter(trend=self, date__lt=from_date)

    def clear_all(self):
        """Delete all matched matches"""
        StreamMatches.objects.filter(stream=self).all().delete()

    def update(self):

        self.set_state('ING')

        # get out-dated StreamMatches
        del_objs = self.get_old_trendmatches()

        # get paper_id list for StreamMatches that remain
        current_paper_list = self.trendmatches_set\
            .exclude(id__in=del_objs)\
            .values_list('paper_id', flat=True)

        # Instantiate Score
        Score = eval(dict(TREND_METHODS_MAP)[self.user.settings.trend_method])
        # and instantiate
        journal_ratio = MostSimilar.objects.get(is_active=True).journal_ratio
        # method_arg = self.user.settings.trend_method_args or {}
        method_arg = {
            'doc_weight': self.user.settings.trend_doc_weight,
            'altmetric_weight': self.user.settings.trend_altmetric_weight,
        }
        scoring = Score(stream=self.user.streams.first(),
                        journal_ratio=journal_ratio,
                        **method_arg)

        # Score
        results, date = scoring.score()

        # create/update UserFeedPaper
        save_objs = []
        create_objs = []
        for pk, val in results:
            if pk in current_paper_list:
                tm = TrendMatches.objects.get(
                    trend=self,
                    paper_id=pk)
                tm.score = val,
                save_objs.append(tm)
            else:
                create_objs.append(TrendMatches(
                    trend=self,
                    paper_id=pk,
                    score=val,
                    date=date[pk]))

        # Get last user visit
        try:
            last_seen = LastSeen.objects.when(user=self.user)
        except LastSeen.DoesNotExist:
            last_seen = None
            pass

        # atomic update
        self.atomic_update(del_objs, save_objs, create_objs, last_seen=last_seen)

        self.set_state('IDL')

    @atomic
    def atomic_update(self, delete_objs, save_objs, create_objs, last_seen=None):
        if delete_objs:
            delete_objs.delete()
        if save_objs:
            for obj in save_objs:
                obj.save()
        if create_objs:
            TrendMatches.objects.bulk_create(create_objs)
        if last_seen:
            TrendMatches.objects.filter(created__lt=last_seen)\
                .update(new=False)


class TrendMatches(TimeStampedModel):

    trend = models.ForeignKey(Trend)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.0)

    date = models.DateField()

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('trend', 'paper'), ]

