# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db import models

from django.contrib.auth.models import BaseUserManager
from django.db import transaction
from django.conf import settings
from django.utils import timezone

from etalia.core.models import TimeStampedModel
from etalia.last_seen.models import LastSeen
from .constants import FEED_STATUS_CHOICES
from etalia.nlp.models import PaperEngine
from etalia.nlp.models import ThreadEngine
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

    def update_async(self):
        from .tasks import update_stream
        update_stream.delay(self.user.id, name=self.name)

    def reset_async(self):
        from .tasks import reset_stream
        reset_stream.delay(self.user.id, name=self.name)

    def update(self):
        from etalia.nlp.tasks import pe_dispatcher

        """Update Stream"""

        logger.info('Updating stream {id}'.format(id=self.id))
        self.set_state('ING')

        # update score threshold
        self.score_threshold = self.user.settings.stream_score_threshold

        # Score
        task = pe_dispatcher.delay('score_stream', self.user.id)
        res = task.get()
        if res:     # res can be empty is user library is empty
            # reformat
            res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                            for r in res if r['score'] > self.score_threshold])
            pids = list(res_dic.keys())

            # clean stream
            self.clean_not_in(pids)

            # Update existing StreamPapers
            with transaction.atomic():
                sp_update = StreamPapers.objects.select_for_update()\
                    .filter(stream=self, paper_id__in=pids)
                update_pids = []
                for sp in sp_update:
                    sp.score = res_dic[sp.paper_id]['score']
                    sp.save()
                    update_pids.append(sp.paper_id)

            # Create new StreamPapers
            create_objs = []
            create_pids = set(pids).difference(set(update_pids))
            for id_ in create_pids:
                create_objs.append(StreamPapers(
                    stream=self,
                    paper_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            StreamPapers.objects.bulk_create(create_objs)

            # Get last user visit. use for the tagging matched with 'new'
            try:
                last_seen = LastSeen.objects.when(user=self.user)
                StreamPapers.objects.filter(created__lt=last_seen).update(new=False)
            except LastSeen.DoesNotExist:
                pass

        self.updated_at = timezone.now()
        self.save(update_fields=('updated_at', ))

        self.set_state('IDL')

        logger.info('Updating stream {id} done'.format(id=self.id))


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

    def update_async(self):
        from .tasks import update_trend
        update_trend.delay(self.user.id, name=self.name)

    def reset_async(self):
        from .tasks import reset_trend
        reset_trend.delay(self.user.id, name=self.name)

    def update(self):
        from etalia.nlp.tasks import pe_dispatcher

        logger.info('Updating trend {id}'.format(id=self.id))
        self.set_state('ING')

        # update score threshold
        self.score_threshold = self.user.settings.stream_score_threshold

        # Score
        task = pe_dispatcher.delay('score_trend', self.user.id)
        res = task.get()
        if res:     # res can be empty is user library is empty
            # reformat
            res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                            for r in res if r['score'] > self.score_threshold])
            pids = list(res_dic.keys())

            # clean trend
            self.clean_not_in(pids)

            # Update existing TrendPapers
            with transaction.atomic():
                tp_update = TrendPapers.objects.select_for_update()\
                    .filter(trend=self, paper_id__in=pids)
                update_pids = []
                for tp in tp_update:
                    tp.score = res_dic[tp.paper_id]['score']
                    tp.save()
                    update_pids.append(tp.paper_id)

            # Create new TrendPapers
            create_objs = []
            create_pids = set(pids).difference(set(update_pids))
            for id_ in create_pids:
                create_objs.append(TrendPapers(
                    trend=self,
                    paper_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            TrendPapers.objects.bulk_create(create_objs)

            # Get last user visit. use for the tagging matched with 'new'
            try:
                last_seen = LastSeen.objects.when(user=self.user)
                TrendPapers.objects.filter(created__lt=last_seen).update(new=False)
            except LastSeen.DoesNotExist:
                pass

        self.updated_at = timezone.now()
        self.save(update_fields=('updated_at', ))

        self.set_state('IDL')
        logger.info('Updating trend {id} done'.format(id=self.id))


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

    def update_async(self):
        from .tasks import update_threadfeed
        update_threadfeed.delay(self.user.id, name=self.name)

    def reset_async(self):
        from .tasks import reset_threadfeed
        reset_threadfeed.delay(self.user.id, name=self.name)

    def update(self):
        from etalia.nlp.tasks import te_dispatcher
        logger.info('Updating thread feed {id}'.format(id=self.id))
        self.set_state('ING')

        # update score threshold
        self.score_threshold = self.user.settings.trend_score_threshold

        # Score
        task = te_dispatcher.delay('score_threadfeed', self.user.id)
        res = task.get()
        if res:     # res can be empty is user library is empty
            # reformat
            res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                            for r in res if r['score'] > self.score_threshold])
            tids = list(res_dic.keys())

            # clean threadfeed
            self.clean_not_in(tids)

            # Update existing TrendFeedThreads
            with transaction.atomic():
                tfp_update = ThreadFeedThreads.objects.select_for_update()\
                    .filter(threadfeed=self, thread_id__in=tids)
                update_tids = []
                for tfp in tfp_update:
                    tfp.score = res_dic[tfp.thread_id]['score']
                    tfp.save()
                    update_tids.append(tfp.thread_id)

            # Create new TrendFeedThreads
            create_objs = []
            create_pids = set(tids).difference(set(update_tids))
            for id_ in create_pids:
                create_objs.append(ThreadFeedThreads(
                    threadfeed=self,
                    thread_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            ThreadFeedThreads.objects.bulk_create(create_objs)

            # Get last user visit. use for the tagging matched with 'new'
            try:
                last_seen = LastSeen.objects.when(user=self.user)
                ThreadFeedThreads.objects.filter(created__lt=last_seen).update(new=False)
            except LastSeen.DoesNotExist:
                pass

        self.updated_at = timezone.now()
        self.save(update_fields=('updated_at', ))

        self.set_state('IDL')
        logger.info('Updating thread feed {id} done'.format(id=self.id))


class ThreadFeedThreads(TimeStampedModel):

    threadfeed = models.ForeignKey(ThreadFeed)

    thread = models.ForeignKey(Thread)

    score = models.FloatField(default=0.0, db_index=True)

    date = models.DateField(db_index=True)

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('threadfeed', 'thread'), ]
