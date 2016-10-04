# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.db import models

from django.conf import settings
from django.core.urlresolvers import reverse
from django.utils import timezone
from django.utils.text import slugify
from django.contrib.auth import get_user_model

from etalia.core.models import TimeStampedModel
from etalia.library.models import Paper
from etalia.core.mixins import ModelDiffMixin
from .constants import THREAD_TYPES, THREADFEED_STATUS_CHOICES, \
    THREAD_TIME_LAPSE_CHOICES, THREAD_INVITE_STATUSES, THREAD_INVITE_PENDING, \
    THREAD_QUESTION, THIRD_PARTY_TYPES, \
    THREAD_PRIVACIES, THREAD_PUBLIC, THREAD_PARTICIPATE, THREAD_WATCH, \
    THREAD_JOINED, THREAD_BANNED, THREAD_PINNED, THREAD_LEFT


class Thread(TimeStampedModel):

    THIRD_PARTY_TYPES = THIRD_PARTY_TYPES

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
    title = models.CharField(max_length=512, verbose_name='Title', default='')

    # content of the thread
    content = models.TextField(null=True, blank=True, default='',
                               verbose_name='Content')

    # published at
    published_at = models.DateTimeField(null=True, blank=True)

    @property
    def short_title(self):
        return self.title[:30]

    @property
    def members(self):
        return self.get_active_members()

    def __str__(self):
        return '{0} | {1} | {2}'.format(self.id,
                                        self.short_title,
                                        self.user)

    def get_absolute_url(self):
        return reverse('threads:thread-slug',
                       kwargs={'pk': self.pk,
                               'slug': slugify(self.title)})

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
        from .tasks import embed_thread
        embed_thread(self.id)

    def get_neighbors(self, time_span):
        from .tasks import get_neighbors_threads
        return get_neighbors_threads(self.id, time_span)


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

    class Meta:
        unique_together = (('thread', 'from_user', 'to_user'), )


class ThreadUser(ModelDiffMixin, TimeStampedModel):
    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # thread
    thread = models.ForeignKey(Thread)

    # Pinned or banned
    watch = models.PositiveIntegerField(null=True, default=None,
                                        choices=THREAD_WATCH)

    # Joined or left
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
        if self.id:
            ThreadUserHistory.objects.create(threaduser=self,
                                             difference=json.dumps(self.diff))
        super(ThreadUser, self).save(**kwargs)


class ThreadUserHistory(TimeStampedModel):

    threaduser = models.ForeignKey(ThreadUser)

    difference = models.CharField(max_length=256, default='')

    date = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ('-created', )


class PubPeer(TimeStampedModel):

    thread = models.OneToOneField(Thread)

    # specify if pubpeer comment should be active in etalia
    is_active = models.BooleanField(default=True)

    doi = models.CharField(max_length=64, blank=True, default='',
                           null=False, unique=True, verbose_name='DOI',
                           db_index=True)

    link = models.URLField()

    pubpeer_id = models.CharField(max_length=30, unique=True)

    @property
    def comments_count(self):
        return self.comments.all().count()

    @property
    def members_count(self):
        members = set([c.user for c in self.comments.all()])
        return len(members)

    @property
    def url(self):
        return self.link

    def update_is_active(self):
        """Activate or deactivate pubpeer thread depending on what is in there

        Currently, deactivating all comments created by bots.
        """
        c = self.comments.first()
        if c.body.startswith('Using the R package statcheck'):
            self.is_active = False
        else:
            self.is_active = True
        self.save(update_fields=['is_active'])


class PubPeerComment(TimeStampedModel):

    pubpeer = models.ForeignKey(PubPeer, related_name='comments')

    body = models.TextField(default='')

    date = models.DateTimeField()

    pubpeercomment_id = models.IntegerField(unique=True)

    permalink = models.URLField()

    rating = models.FloatField()

    user = models.CharField(max_length=128)