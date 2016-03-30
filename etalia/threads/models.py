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
    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='threads_owned')

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

    class Meta:
        unique_together = (('type', 'user', 'title', 'paper'), )

    @property
    def short_title(self):
        return self.title[:30]

    @property
    def members(self):
        return self.get_active_members()

    def __str__(self):
        return '{0}@{1}'.format(self.short_title, self.user)

    def get_active_members(self):
        tus = self.threaduser_set.filter(is_joined=True).select_related('user')
        return [tu.user for tu in tus]

    def is_owner(self, user):
        if user == self.user:
            return True
        else:
            return False

    def save(self, *args, **kwargs):
        if not self.threaduser_set.exists():
            tu = self.threaduser_set.create(user=self.user)
            tu.join()
        super(Thread, self).save(*args, **kwargs)

    @property
    def state(self, user):
        if ThreadUser.objects.filter(user=user, thread=self).exists():
            return ThreadUser.objects.get(user=user, thread=self)
        else:
            return None


class ThreadPost(TimeStampedModel):

    # thread
    thread = models.ForeignKey(Thread, related_name='posts')

    # User who posts
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('thread', 'content', 'user'), )
        ordering = ('created', )

    def __str__(self):
        return '{0} {1}'.format(self.created, self.user)

    def save(self, **kwargs):
        if not self.id:
            # increment ThreadUser post count
            self.user.threaduser_set\
                .filter(thread=self.thread)\
                .update(num_comments=F('num_comments') + 1)
        super(ThreadPost, self).save(**kwargs)

    def delete(self, **kwargs):
        # decrement ThreadUser post count
        self.user.threaduser_set\
            .filter(thread=self.thread)\
            .update(num_comments=F('num_comments') - 1)
        super(ThreadPost, self).delete(**kwargs)


class ThreadComment(TimeStampedModel):

    # thread
    post = models.ForeignKey(ThreadPost, related_name='comments')

    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # content
    content = models.TextField(null=False, blank=True, default='')

    class Meta:
        unique_together = (('post', 'user', 'content'), )
        ordering = ('created', )

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


class ThreadUser(TimeStampedModel):

    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # thread
    thread = models.ForeignKey(Thread)

    # thread is pinned
    is_pinned = models.BooleanField(default=False)

    # thread is banned
    is_banned = models.BooleanField(default=False)

    # thread is added
    is_joined = models.BooleanField(default=False)

    # thread is trashed
    is_left = models.BooleanField(default=False)

    # First time joined the thread
    first_joined_at = models.DateTimeField(null=True, blank=True, default=None)

    # When user left the thread, is any
    last_left_at = models.DateTimeField(null=True, blank=True, default=None)

    # Number of comments posted in thread
    num_comments = models.PositiveIntegerField(default=0)

    class Meta:
        ordering = ['-num_comments', 'first_joined_at']

    def join(self):
        self.is_joined = True
        self.is_left = False
        self.is_banned = False
        if not self.first_joined_at:
            self.first_joined_at = timezone.now()
        self.last_left_at = None
        self.save()

    def leave(self):
        self.is_joined = False
        self.is_banned = False
        self.is_left = True
        self.last_left_at = timezone.now()
        self.save()

    def pin(self):
        self.is_pinned = not self.is_pinned
        self.is_banned = False
        self.save()

    def ban(self):
        self.is_pinned = False
        self.is_banned = True
        self.save()



