# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from celery.canvas import chain

from django.db import models


from django.contrib.auth.models import BaseUserManager
from django.db import transaction
from django.conf import settings
from django.utils import timezone

from etalia.core.models import TimeStampedModel
from etalia.last_seen.models import LastSeen
from .constants import FEED_STATUS_CHOICES
from etalia.threads.models import Thread
from etalia.library.models import Paper

logger = logging.getLogger(__name__)


class StreamManager(BaseUserManager):

    def create(self, **kwargs):
        """Create new UserFeed / UserFeedVector
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')

        obj = super(StreamManager, self).create(**kwargs)

        return obj

    def create_main(self, **kwargs):
        """Create a new default UserFeed / UserFeedVector, named 'main' """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')
        return self.create(name='main', **kwargs)


class Stream(TimeStampedModel):
    """Stream: Table of relevant Papers for a user"""

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='streams')

    name = models.CharField(max_length=128, default='main')

    papers = models.ManyToManyField(Paper, through='StreamPapers',
                                    related_name='stream_papers')

    state = models.CharField(max_length=3, blank=True, default='NON',
                             choices=FEED_STATUS_CHOICES)

    score_threshold = models.FloatField(default=0.3)

    updated_at = models.DateTimeField(default=None, blank=True, null=True)

    objects = StreamManager()

    class Meta:
        unique_together = (('name', 'user'),)

    def set_state(self, state):
        self.state = state
        self.save(update_fields=('state', ))

    def __str__(self):
        return '{stream}@{email}'.format(stream=self.name,
                                         email=self.user.email)

    def clear_all(self):
        """Delete all papers"""
        StreamPapers.objects.filter(stream=self).delete()

    def clean_not_in(self, paper_ids):
        """Delete papers that are not in paper_ids"""
        StreamPapers.objects\
            .filter(stream=self)\
            .exclude(paper_id__in=paper_ids)\
            .delete()

    def update(self):
        """Update Stream"""
        from etalia.nlp.tasks import pe_dispatcher
        from .tasks import populate_stream

        # Chain tasks
        chain(pe_dispatcher.s('score_stream', self.id),
              populate_stream.s(self.id))()


class StreamPapers(TimeStampedModel):
    """Relationship table between stream and matched matches with scoring"""

    stream = models.ForeignKey(Stream)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.0, db_index=True)

    date = models.DateField(db_index=True)

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('stream', 'paper')]

    def print_score(self):
        return '{0:.1f}'.format(self.score * 100)

    def __str__(self):
        return '{paper}/{score}'.format(paper=self.paper.short_title,
                                        score=self.score)


class Trend(TimeStampedModel):
    """Trend: Table of relevant trendy Papers for a user"""

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='trends')

    name = models.CharField(max_length=128, default='main')

    papers = models.ManyToManyField(Paper, through='TrendPapers')

    state = models.CharField(max_length=3, blank=True, choices=FEED_STATUS_CHOICES)

    updated_at = models.DateTimeField(default=None, blank=True, null=True)

    score_threshold = models.FloatField(default=0.1)

    def __str__(self):
        return self.user.email

    class Meta:
        unique_together = (('user', 'name'),)

    def set_state(self, state):
        self.state = state
        self.save(update_fields=('state', ))

    def clear_all(self):
        """Delete all matched matches"""
        TrendPapers.objects.filter(trend=self).all().delete()

    def clean_not_in(self, paper_ids):
        """Delete papers that are not in paper_ids"""
        TrendPapers.objects\
            .filter(trend=self)\
            .exclude(paper_id__in=paper_ids)\
            .delete()

    def update(self):
        """Update Trend"""
        from etalia.nlp.tasks import pe_dispatcher
        from .tasks import populate_trend

        # Chain tasks
        chain(pe_dispatcher.s('score_trend', self.id),
              populate_trend.s(self.id))()


class TrendPapers(TimeStampedModel):

    trend = models.ForeignKey(Trend)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.0, db_index=True)

    date = models.DateField(db_index=True)

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('trend', 'paper'), ]


class ThreadFeed(TimeStampedModel):
    """Trend: Table of relevant Threads for a user"""

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='threadfeeds')

    name = models.CharField(max_length=128, default='main')

    threads = models.ManyToManyField(Thread, through='ThreadFeedThreads')

    state = models.CharField(max_length=3, blank=True, choices=FEED_STATUS_CHOICES)

    updated_at = models.DateTimeField(default=None, blank=True, null=True)

    score_threshold = models.FloatField(default=0.3)

    def __str__(self):
        return self.user.email

    class Meta:
        unique_together = (('user', 'name'),)

    def set_state(self, state):
        self.state = state
        self.save()

    def clear_all(self):
        """Delete all matched matches"""
        ThreadFeedThreads.objects.filter(threadfeed=self).all().delete()

    def clean_not_in(self, thread_ids):
        """Delete threads that are not in thread_ids"""
        ThreadFeedThreads.objects\
            .filter(threadfeed=self)\
            .exclude(thread_id__in=thread_ids)\
            .delete()

    def update(self):
        from etalia.nlp.tasks import te_dispatcher
        from .tasks import populate_threadfeed

        # Chain tasks
        chain(te_dispatcher.s('score_threadfeed', self.id),
              populate_threadfeed.s(self.id))()


class ThreadFeedThreads(TimeStampedModel):

    threadfeed = models.ForeignKey(ThreadFeed)

    thread = models.ForeignKey(Thread)

    score = models.FloatField(default=0.0, db_index=True)

    date = models.DateField(db_index=True)

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('threadfeed', 'thread'), ]
