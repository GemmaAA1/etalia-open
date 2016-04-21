# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.db import models

from django.conf import settings
from django.utils import timezone
from django.contrib.auth import get_user_model
from django.contrib.postgres.fields import ArrayField

from etalia.core.models import TimeStampedModel


from etalia.library.models import Paper
from .mixins import ModelDiffMixin
from .constant import THREAD_TYPES, THREADFEED_STATUS_CHOICES, \
    THREAD_TIME_LAPSE_CHOICES, THREAD_INVITE_STATUSES, THREAD_INVITE_PENDING, \
    THREAD_QUESTION, \
    THREAD_PRIVACIES, THREAD_PUBLIC, THREAD_PARTICIPATE, THREAD_WATCH, \
    THREAD_JOINED, THREAD_BANNED, THREAD_PINNED, THREAD_LEFT


class Thread(TimeStampedModel):

    # type of thread
    type = models.IntegerField(choices=THREAD_TYPES, default=THREAD_QUESTION,
                               null=False, blank=False, verbose_name='Type')

    # User who create thread
    user = models.ForeignKey(settings.AUTH_USER_MODEL,
                             related_name='threads_owned',
                             null=True)

    # Privacy
    privacy = models.IntegerField(choices=THREAD_PRIVACIES,
                                  default=THREAD_PUBLIC)

    # Paper that the thread is based on, if any
    paper = models.ForeignKey(Paper, null=True, blank=True, default=None,
                              verbose_name='Related Paper')

    # title of thread
    title = models.CharField(max_length=256, verbose_name='Title', default='')

    # content of the thread
    content = models.TextField(null=True, blank=True, default='',
                               verbose_name='Content')

    # published at
    published_at = models.DateTimeField(null=True, blank=True)

    class Meta:
        ordering = ('-published_at', '-modified')

    @property
    def short_title(self):
        return self.title[:30]

    @property
    def members(self):
        return self.get_active_members()

    def __str__(self):
        return '{0}@{1}'.format(self.short_title, self.user)

    def get_active_members(self):
        User = get_user_model()
        return User.objects.filter(threaduser__thread_id=self.id,
                                   threaduser__participate=THREAD_JOINED)

    def is_owner(self, user):
        if user == self.user:
            return True
        else:
            return False

    def save(self, *args, **kwargs):
        super(Thread, self).save(*args, **kwargs)
        if not self.threaduser_set.exists():
            tu = self.threaduser_set.create(user=self.user)
            tu.join()

    def state(self, user):
        if ThreadUser.objects.filter(user=user, thread=self).exists():
            return ThreadUser.objects.get(user=user, thread=self)
        else:
            return None

    def publish(self):
        self.published_at = timezone.now()
        self.save(update_fields=['published_at'])

    def embed(self):
        pass

    def get_members_per_post_count(self):
        pass


class ThreadPost(TimeStampedModel):
    # thread
    thread = models.ForeignKey(Thread, related_name='posts')

    # User who posts
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('thread', 'content', 'user'),)
        ordering = ('created',)

    def __str__(self):
        return '{0} {1}'.format(self.created, self.user)


class ThreadComment(TimeStampedModel):
    # thread
    post = models.ForeignKey(ThreadPost, related_name='comments')

    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('post', 'user', 'content'),)
        ordering = ('created',)

    def __str__(self):
        return '{0} {1}'.format(self.created.ctime(), self.user)


class ThreadUserInvite(TimeStampedModel):
    # thread
    thread = models.ForeignKey(Thread)

    # user
    from_user = models.ForeignKey(settings.AUTH_USER_MODEL,
                                  related_name='invite_from_users')

    # invite user
    to_user = models.ForeignKey(settings.AUTH_USER_MODEL,
                                related_name='invite_to_users')

    # invite status
    status = models.IntegerField(choices=THREAD_INVITE_STATUSES,
                                 default=THREAD_INVITE_PENDING)


class ThreadFeed(TimeStampedModel):
    # User
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # name of the feed
    name = models.CharField(max_length=128, default='main')

    # last update of the feed
    updated_at = models.DateTimeField(auto_now=True)

    # state
    status = models.CharField(max_length=3, blank=True,
                              default='NON',
                              choices=THREADFEED_STATUS_CHOICES)

    # threads
    thread_scores = models.ManyToManyField(Thread, through='ThreadScore')


class ThreadScore(TimeStampedModel):
    # thread
    thread = models.ForeignKey(Thread)

    # feed
    thread_feed = models.ForeignKey(ThreadFeed)

    # score
    score = models.FloatField(default=0.0)

    # new flag
    new = models.BooleanField(default=True)

    class Meta:
        ordering = ['-score']
        unique_together = (('thread', 'thread_feed'),)


class ThreadUser(ModelDiffMixin, TimeStampedModel):
    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # thread
    thread = models.ForeignKey(Thread)

    # thread watch (pinned or left
    watch = models.PositiveIntegerField(null=True, default=None,
                                        choices=THREAD_WATCH)

    # thread is banned
    participate = models.PositiveIntegerField(null=True, default=None,
                                              choices=THREAD_PARTICIPATE)

    class Meta:
        unique_together = (('thread', 'user'),)

    def join(self):
        self.participate = THREAD_JOINED
        self.save()

    def leave(self):
        self.participate = THREAD_LEFT
        self.save()

    def pin(self):
        self.watch = THREAD_PINNED
        self.save()

    def ban(self):
        self.watch = THREAD_BANNED
        self.save()

    def save(self, **kwargs):
        super(ThreadUser, self).save(**kwargs)
        ThreadUserHistory.objects.create(threaduser=self,
                                         difference=json.dumps(self.diff))


class ThreadUserHistory(TimeStampedModel):
    threaduser = models.ForeignKey(ThreadUser)

    difference = models.CharField(max_length=256, default='')

    date = models.DateTimeField(auto_now_add=True)


class ThreadNeighbor(TimeStampedModel):
    """Table for Thread neighbors"""
    # thread
    thread = models.ForeignKey(Thread)

    # Time lapse for neighbors match
    time_lapse = models.IntegerField(default=-1,
                                     choices=THREAD_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)
