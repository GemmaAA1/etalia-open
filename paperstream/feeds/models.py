# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

import numpy as np
from django.db import models
from django.conf import settings
from django.contrib.auth.models import BaseUserManager
from django.db.models import Q
from django.utils import timezone
from django.contrib.postgres.fields import ArrayField

from gensim import matutils

from paperstream.core.models import TimeStampedModel
from paperstream.core.utils import pad_vector
from paperstream.library.models import Paper
from paperstream.nlp.models import Model, MostSimilar
from config.celery import celery_app as app
from .constants import FEED_STATUS_CHOICES
from .utils import SimpleAverage, ThresholdAverage, WeightedJournalAverage, \
    WeightedJournalCreatedDateAverage, SimpleMax

logger = logging.getLogger(__name__)


class UserFeedManager(BaseUserManager):

    def create(self, **kwargs):
        """Create new UserFeed / UserFeedVector
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')

        papers_seed = None
        if 'papers_seed' in kwargs:
            papers_seed = kwargs.pop('papers_seed')

        obj = super(UserFeedManager, self).create(**kwargs)

        if papers_seed:
            obj.add_papers_seed(papers_seed)
            obj.save()

        return obj

    def create_main(self, **kwargs):
        """Create a new default UserFeed / UserFeedVector, named 'main' """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')
        return self.create(name='main', **kwargs)


class UserFeed(TimeStampedModel):
    """User Feed model"""

    # feed name
    name = models.CharField(max_length=100, default='main')

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='feeds')

    # Papers that are used in the similarity matching
    papers_seed = models.ManyToManyField(Paper, through='UserFeedSeedPaper',
                                         related_name='feed_seed')

    # relevant papers matched
    papers_match = models.ManyToManyField(Paper, through='UserFeedMatchPaper',
                                          related_name='feed_match')

    state = models.CharField(max_length=3, blank=True, default='NON',
                             choices=FEED_STATUS_CHOICES)

    message = models.CharField(max_length=127, blank=True, default='')

    objects = UserFeedManager()

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
        return '{feed}@{username}'.format(feed=self.name, username=self.user.email)


    def add_papers_seed(self, papers):

        if papers:
            objs = []
            for paper in papers:
                objs.append(UserFeedSeedPaper(feed=self,
                                              paper=paper))
            UserFeedSeedPaper.objects.bulk_create(objs)

        # Update feed vector
        self.update_userfeed_vector()

    def update_userfeed_vector(self):
        model_pks = Model.objects.all().values_list('pk', flat='True')
        for model_pk in model_pks:
            ufv, _ = UserFeedVector.objects.get_or_create(model_id=model_pk,
                                                          feed=self)
            ufv.update_vector()

    def clean_old_papers(self):
        # get time lapse corresponding to user settings
        from_date = (timezone.now() - timezone.timedelta(
            days=self.user.settings.time_lapse)).date()

        # clean old papers
        UserFeedMatchPaper.objects\
            .filter(Q(feed=self) &
                    (Q(paper__date_ep__lt=from_date) |
                    (Q(paper__date_pp__lt=from_date) & Q(paper__date_ep=None))))\
            .delete()

    def clear(self):
        UserFeedMatchPaper.objects.filter(feed=self).delete()

    def log(self, level, core, message):
        if level == 'info':
            logger.info('Feed ({pk}/{feed_name}@{user_email}): {core} - {message}'.format(
                pk=self.id, feed_name=self.name, user_email=self.user.email, core=core, message=message))
        elif level == 'debug':
            logger.debug('Feed ({pk}/{feed_name}@{user_email}): {core} - {message}'.format(
                pk=self.id, feed_name=self.name, user_email=self.user.email, core=core, message=message))
        else:
            raise ValueError('level unknown')

    def update(self, restrict_journal=False):
        """Update UserFeed

        - Get all papers in time range excluding untrusted and already in user lib
        #TODO: Restricted this query to the k-NN to the feed vector
        - Score papers
        - Get only the N top scored papers
        - Create/Update UserFeedPaper
        """


        self.set_state('ING')
        self.log('info', 'Updating', 'starting...')
        self.set_message('Cleaning')

        # self.clean_old_papers()
        self.clear()

        # get paper that remain and need to be updated
        self.set_message('Fetching data')
        ufp_pks_to_update = self.papers_match.values_list('pk', flat=True)
        # get user lib pks
        lib_pks = self.user.lib.papers.values_list('pk', flat=True)

        # instantiate scoring
        if self.user.settings.scoring_method == 0:
            Score = SimpleMax
        elif self.user.settings.scoring_method == 1:
            Score = SimpleAverage
        elif self.user.settings.scoring_method == 2:
            Score = ThresholdAverage
        elif self.user.settings.scoring_method == 3:
            Score = WeightedJournalAverage
        elif self.user.settings.scoring_method == 4:
            Score = WeightedJournalCreatedDateAverage
        else:
            raise ValueError

        scoring = Score(
            model=self.user.settings.model,
            user=self.user)

        # get target papers (excluding userlib + not trusted + empty abstract)
        # from MostSimilar
        # Get corresponding task:
        try:
            ms_task = app.tasks['paperstream.nlp.tasks.mostsimilar_{name}'.format(
                name=self.user.settings.model.name)]
        except KeyError:
            raise KeyError

        self.log('debug', 'Updating', 'getting neighbors...')
        seed_pks = self.papers_seed.all().values_list('pk', flat=True)
        # Submit Celery Task to get k_neighbors
        res = ms_task.delay('get_partition',
                             paper_pks=seed_pks,
                             time_lapse=self.user.settings.time_lapse,
                             k=settings.FEED_K_NEIGHBORS)
        # Wait for Results
        target_seed_pks = res.get()
        self.log('debug', 'Updating', 'done')
        # Filter exclude corresponding papers
        target_pks = list(ufp_pks_to_update) + \
                     [pk for pk in target_seed_pks if pk not in lib_pks]

        # Match target paper
        objs_list = []
        if target_pks:
            # compute scores
            self.log('debug', 'Updating', 'preparing...')
            scoring.prepare(seed_pks, target_pks)
            self.log('debug', 'Updating', 'done')
            self.set_message('Scoring {0} papers'.format(len(target_pks)))
            self.log('debug', 'Updating', 'scoring...')
            pks, scores = scoring.score()
            self.log('debug', 'Updating', 'done')
            # sort scores
            self.log('debug', 'Updating', 'sorting...')
            best = matutils.argsort(scores,
                                    topn=settings.FEEDS_SCORE_KEEP_N_PAPERS,
                                    reverse=True)
            # reshape
            results = [(pks[ind], float(scores[ind])) for ind in best]
            self.log('debug', 'Updating', 'done')

            # create/update UserFeedPaper
            self.log('debug', 'Updating', 'populating...')
            for pk, val in results:
                # update
                if pk in ufp_pks_to_update:
                    ufp, _ = UserFeedMatchPaper.objects.get_or_create(
                        feed=self,
                        paper_id=pk)
                    ufp.score = val
                    ufp.is_score_computed = True
                    ufp.save()
                else:  # create in bulk
                    objs_list.append(UserFeedMatchPaper(
                        feed=self,
                        paper_id=pk,
                        score=val,
                        is_score_computed=True))
            # bulk create
            UserFeedMatchPaper.objects.bulk_create(objs_list)
            self.log('debug', 'Updating', 'done')
        self.set_state('IDL')

        self.log('info', 'Updating', 'DONE')


class UserFeedSeedPaper(TimeStampedModel):
    """Relationship table between feed and seed papers"""

    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    class Meta:
        unique_together = ('feed', 'paper')

    def clean(self):
        """check if paper is in user.lib"""
        assert self.paper.pk in \
            self.feed.user.lib.papers.values_list('pk', flat=True)


class UserFeedMatchPaper(TimeStampedModel):
    """Relationship table between feed and matched papers with scoring"""
    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    is_score_computed = models.BooleanField(default=False)

    class Meta:
        unique_together = ('feed', 'paper')

    class Meta:
        ordering = ['-score']
        unique_together = [('feed', 'paper')]

    def print_score(self):
        return '{0:.1f}'.format(self.score * 100)

    def __str__(self):
        return '{paper}/{score}'.format(paper=self.paper.short_title,
                                        score=self.score)

    @property
    def is_disliked(self):
        return UserTaste.objects\
            .get(paper=self.paper, user=self.feed.user).is_disliked

    @property
    def is_liked(self):
        return UserTaste.objects\
            .get(paper=self.paper, user=self.feed.user).is_liked

class UserFeedVector(TimeStampedModel):
    """Feature vector for feed is defined as the averaged of paper vectors in
    feed
    """
    feed = models.ForeignKey(UserFeed, related_name='vectors')

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True)

    class Meta:
        unique_together = ('feed', 'model')

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        if self.vector:
            return self.vector[:self.model.size]
        else:
            return None

    def update_vector(self):
        vectors = self.feed.papers_seed\
            .filter(vectors__model=self.model)\
            .values_list('vectors__vector', flat=True)

        # sum of vectors
        vector = [sum(x) for x in zip(*vectors) if x[0]]  # not considering None
        self.set_vector(vector)

    class Meta:
        unique_together = [('feed', 'model'), ]

    def __str__(self):
        return '{feed_name}/{model_name}'.format(feed_name=self.feed.name,
                                                 model_name=self.model.name)

