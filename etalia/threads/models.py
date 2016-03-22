# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.conf import settings

from etalia.core.models import TimeStampedModel
from etalia.library.models import Paper
from .constant import THREAD_TYPES, THREADFEED_STATUS_CHOICES, \
    THREAD_TIME_LAPSE_CHOICES, INVITE_STATUSES, INVITE_PENDING, THREAD_QUESTION


class Thread(TimeStampedModel):

    # type of thread
    type = models.IntegerField(choices=THREAD_TYPES, default=THREAD_QUESTION,
                               null=False, blank=False)

    # User who create thread
    author = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='threads')

    # Users that joined the thread
    members = models.ManyToManyField(settings.AUTH_USER_MODEL,
                                     through='ThreadMember')

    # Paper that the thread is based on, if any
    paper = models.ForeignKey(Paper, null=True, blank=True, default=None)

    # title of thread
    title = models.CharField(max_length=256)

    # content of the thread
    content = models.TextField(null=True, blank=True, default='')


class ThreadMember(models.Model):

    # thread
    thread = models.ForeignKey(Thread)

    # member of the thread
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # First time joined the thread
    joined_at = models.DateTimeField(auto_now_add=True)

    # When user left the thread, is any
    left_at = models.DateTimeField(null=True, blank=True, default=None)

    # Number of comments posted in thread
    num_comments = models.PositiveIntegerField(default=0)


class ThreadPost(TimeStampedModel):

    # thread
    thread = models.ForeignKey(Thread, related_name='posts')

    # User who posts
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # index of the post in the thread
    position = models.PositiveIntegerField(default=0)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('thread', 'position'), )


class ThreadPostComment(TimeStampedModel):

    # thread
    post = models.ForeignKey(ThreadPost, related_name='comments')

    # index of the comment for this post
    position = models.PositiveIntegerField(default=0)

    # author
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('post', 'position'), )


class ThreadUserState(TimeStampedModel):

    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # thread
    thread = models.ForeignKey(Thread)

    # thread is pinned
    is_pinned = models.BooleanField(default=False)

    # thread is banned
    is_banned = models.BooleanField(default=False)

    # thread is added
    is_added = models.BooleanField(default=False)

    # thread is trashed
    is_trashed = models.BooleanField(default=False)


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
    status = models.IntegerField(choices=INVITE_STATUSES,
                                 default=INVITE_PENDING)


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
    threads = models.ManyToManyField(Thread, through='ThreadFeedThread')


class ThreadFeedThread(TimeStampedModel):

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
        unique_together = (('thread', 'thread_feed'), )


class ThreadNeighbor(TimeStampedModel):

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



