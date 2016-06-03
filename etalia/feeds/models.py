# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db import models

from django.contrib.auth.models import BaseUserManager
from django.db import transaction

from etalia.core.models import TimeStampedModel
from etalia.last_seen.models import LastSeen
from .constants import FEED_STATUS_CHOICES, STREAM_METHODS_MAP, TREND_METHODS_MAP
from .scoring import *
from etalia.nlp.models import PaperEngine


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

    last_update = models.DateTimeField(default=None, blank=True, null=True)

    objects = StreamManager()

    class Meta:
        unique_together = (('name', 'user'),)

    def set_state(self, state):
        self.state = state
        self.save()

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

        logger.info('Updating stream {id}'.format(id=self.id))
        self.set_state('ING')

        # Score
        pe = PaperEngine.objects.get(is_active=True)
        pe_tasks = app.tasks['etalia.nlp.tasks.{name}'.format(name=pe.name)]
        task = pe_tasks.delay('score_stream', self.user.id)
        res = task.get()
        # reformat
        res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']}) for r in res])
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

        self.last_update = timezone.now()
        self.save(update_fields=('last_update', ))

        self.set_state('IDL')
        logger.info('Updating stream {id} done'.format(id=self.id))


class StreamPapers(TimeStampedModel):
    """Relationship table between stream and matched matches with scoring"""

    stream = models.ForeignKey(Stream)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    date = models.DateField()

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

    last_update = models.DateTimeField(default=None, blank=True, null=True)

    def __str__(self):
        return self.user.email

    class Meta:
        unique_together = (('user', 'name'),)

    def set_state(self, state):
        self.state = state
        self.save()

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

        logger.info('Updating trend {id}'.format(id=self.id))
        self.set_state('ING')

        # Score
        pe = PaperEngine.objects.get(is_active=True)
        pe_tasks = app.tasks['etalia.nlp.tasks.{name}'.format(name=pe.name)]
        task = pe_tasks.delay('score_trend', self.user.id)
        res = task.get()
        # reformat
        res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']}) for r in res])
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

        self.last_update = timezone.now()
        self.save(update_fields=('last_update', ))

        self.set_state('IDL')
        logger.info('Updating trend {id} done'.format(id=self.id))


class TrendPapers(TimeStampedModel):

    trend = models.ForeignKey(Trend)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.0)

    date = models.DateField()

    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = [('trend', 'paper'), ]

