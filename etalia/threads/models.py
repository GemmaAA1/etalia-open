# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from django.db.models import F
from django.utils import timezone

from etalia.core.models import TimeStampedModel
from etalia.library.models import Paper
from .constant import THREAD_TYPES, THREADFEED_STATUS_CHOICES, \
    THREAD_TIME_LAPSE_CHOICES, THREAD_INVITE_STATUSES, THREAD_INVITE_PENDING, THREAD_QUESTION,\
    THREAD_PRIVACIES, THREAD_PUBLIC


class Thread(TimeStampedModel):

    # type of thread
    type = models.IntegerField(choices=THREAD_TYPES, default=THREAD_QUESTION,
                               null=False, blank=False, verbose_name='Type')

    # User who create thread
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='threads_owned')

    # Privacy
    privacy = models.IntegerField(choices=THREAD_PRIVACIES,
                                  default=THREAD_PUBLIC)

    # Paper that the thread is based on, if any
    paper = models.ForeignKey(Paper, null=True, blank=True, default=None,
                              verbose_name='Related Paper')

    # title of thread
    title = models.CharField(max_length=256, verbose_name='Title')

    # content of the thread
    content = models.TextField(null=True, blank=True, default='',
                               verbose_name='Content')

    class Meta:
        unique_together = (('type', 'owner', 'title', 'paper'), )

    @property
    def short_title(self):
        return self.title[:30]

    @property
    def members(self):
        return list(self.get_active_members())

    def __str__(self):
        return '{0}@{1}'.format(self.short_title, self.owner)

    def get_active_members(self):
        return self.user_set.filter(userthread__is_joined=True)

    def is_owner(self, user):
        if user == self.owner:
            return True
        else:
            return False

    def save(self, *args, **kwargs):
        super(Thread, self).save(*args, **kwargs)
        # add owner as member
        self.owner.join_thread(self)


class ThreadPost(TimeStampedModel):

    # thread
    thread = models.ForeignKey(Thread, related_name='posts')

    # User who posts
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('thread', 'content', 'author'), )
        ordering = ('created', )

    def __str__(self):
        return '{0} {1}'.format(self.created, self.author)

    def save(self, **kwargs):
        if not self.id:
            # increment UserThread post count
            self.author.userthread_set\
                .filter(thread=self.thread)\
                .update(num_comments=F('num_comments') + 1)
        super(ThreadPost, self).save(**kwargs)

    def delete(self, **kwargs):
        # decrement UserThread post count
        self.author.userthread_set\
            .filter(thread=self.thread)\
            .update(num_comments=F('num_comments') - 1)
        super(ThreadPost, self).delete(**kwargs)


class ThreadPostComment(TimeStampedModel):

    # thread
    post = models.ForeignKey(ThreadPost, related_name='comments')

    # user
    author = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('post', 'author', 'content'), )
        ordering = ('created', )

    def __str__(self):
        return '{0} {1}'.format(self.created.ctime(), self.author)


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



