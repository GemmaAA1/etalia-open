# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
from collections import Counter

from django.db import models
from django.db.models import F, Min
from django.utils import timezone
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.users.models import UserLibPaper
from etalia.library.models import AuthorPaper
from etalia.threads.models import ThreadUser, Thread
from etalia.threads.constant import THREAD_JOINED
from etalia.core.models import TimeStampedModel
from etalia.core.utils import pad_or_trim_vector

from ..constants import FINGERPRINT_STATUS_CHOICES


class UserFingerprintManager(models.Manager):

    def create(self, **kwargs):
        obj = super(UserFingerprintManager, self).create(**kwargs)

        if not obj.user.settings.fingerprint_roll_back_deltatime and \
                        obj.user.lib.papers.count() > 0:
                obj.user.settings.init_fingerprint_roll_back_deltatime()
                obj.update_added_after()

        return obj


class UserFingerprint(TimeStampedModel):

    user = models.ForeignKey(settings.AUTH_USER_MODEL,
                             related_name='fingerprint')

    name = models.CharField(max_length=128, default='main')

    model = models.ForeignKey('nlp.Model')

    added_after = models.DateField(null=True, blank=True)

    embedding = ArrayField(models.FloatField(null=True),
                           size=settings.NLP_MAX_VECTOR_SIZE,
                           null=True)

    authors_ids = ArrayField(models.PositiveIntegerField(null=True),
                             size=settings.NLP_MAX_VECTOR_SIZE,
                             null=True)

    authors_counts = ArrayField(models.PositiveIntegerField(null=True),
                                size=settings.NLP_MAX_VECTOR_SIZE,
                                null=True)

    journals_ids = ArrayField(models.PositiveIntegerField(null=True),
                              size=settings.NLP_MAX_VECTOR_SIZE,
                              null=True)

    journals_counts = ArrayField(models.PositiveIntegerField(null=True),
                                 size=settings.NLP_MAX_VECTOR_SIZE,
                                 null=True)

    thread_embedding = ArrayField(models.FloatField(null=True),
                                  size=settings.NLP_MAX_VECTOR_SIZE,
                                  null=True)

    threads_users_ids = ArrayField(models.PositiveIntegerField(null=True),
                                   size=settings.NLP_MAX_VECTOR_SIZE,
                                   null=True)

    threads_users_counts = ArrayField(models.PositiveIntegerField(null=True),
                                      size=settings.NLP_MAX_VECTOR_SIZE,
                                      null=True)

    state = models.CharField(max_length=3, blank=True,
                             choices=FINGERPRINT_STATUS_CHOICES)

    objects = UserFingerprintManager()

    class Meta:
        unique_together = (('user', 'name'), )

    def __str__(self):
        return '{email}'.format(email=self.user.email)

    @property
    def embedding_size(self):
        return self.model.size

    def set_state(self, state):
        self.state = state
        self.save(update_fields=('state', ))

    def update_added_after(self):
        delta = self.user.settings.fingerprint_roll_back_deltatime  # in months
        added_after = (timezone.now() - timezone.timedelta(days=delta*30.5)).date()
        if not self.added_after == added_after:
            self.added_after = added_after
            self.save(update_fields=('added_after', ))
        return self.added_after

    def update_async(self):
        from ..tasks import update_userfingerprint
        update_userfingerprint.delay(self.user.id, name=self.name)

    def update(self):

        if self.user.lib.papers.count() > 0:

            self.set_state('ING')

            data = {'ids': [],
                    'journal-ids': [],
                    'authors-ids': [],
                    'date': [],
                    'embedding': [],
                    'thread-embedding': [],
                    'thread-users-ids': [],
                    }

            self.update_added_after()

            # Paper related fingerprint
            q1 = UserLibPaper.objects.raw(
                    "SELECT ulp.id, "
                    "		ulp.paper_id, "
                    "		ulp.date_created AS date_, "
                    "		lp.journal_id, "
                    "		pv.vector "
                    "FROM users_userlibpaper ulp "
                    "LEFT JOIN nlp_papervectors pv ON ulp.paper_id = pv.paper_id "
                    "LEFT JOIN library_paper lp ON ulp.paper_id = lp.id "
                    "WHERE ulp.userlib_id = %s "
                    "    AND ulp.date_created >= %s "
                    "    AND pv.vector IS NOT NULL "
                    "ORDER BY ulp.date_created ASC", (self.user_id, self.added_after)
                    )

            for d in q1:
                data['ids'].append(d.paper_id)
                data['journal-ids'].append(d.journal_id)
                data['date'].append(d.date_)
                data['embedding'].append(d.vector[:self.embedding_size])
            data['embedding'] = np.array(data['embedding'])

            # query on authors
            if data['ids']:
                values = ', '.join(['({0})'.format(i) for i in data['ids']])
                q2 = AuthorPaper.objects.raw(
                        "SELECT ap.id, "
                        "	    ap.paper_id, "
                        "		ap.author_id "
                        "FROM library_authorpaper ap "
                        "WHERE ap.paper_id IN (VALUES {0}) ".format(values))

                d2 = [(d.paper_id, d.author_id) for d in q2]
                dic = {}
                for k, v in d2:
                    try:
                        dic[k].append(v)
                    except KeyError:
                        dic[k] = [v]

                for i in data['ids']:
                    if i in dic.keys():
                        data['authors-ids'].append(dic[i])
                    else:
                        data['authors-ids'].append([])

                # Compute date cutoff
                start_ind = np.argmax(np.array(data['date']) > self.added_after)

                # Count (order from higher count to lower count)
                auth_count = Counter(
                    [auth for auths in data['authors-ids'][start_ind:]
                     for auth in auths]).most_common()
                jour_count = Counter(data['journal-ids'][start_ind:]).most_common()

                # Normalize averaged embedding
                vec_paper = np.sum(data['embedding'][start_ind:, :], axis=0)
                norm = np.linalg.norm(vec_paper)
                if norm > 0:
                    vec_paper /= norm

                # Thread related fingerprint
                q2 = Thread.objects.raw(
                    "SELECT t.id, "
                    "		t.user_id,"
                    "		tv.vector "
                    "FROM threads_thread t "
                    "LEFT JOIN nlp_threadvectors tv ON t.id = tv.thread_id "
                    "LEFT JOIN threads_threaduser tu ON t.id = tu.thread_id "
                    "WHERE tu.user_id = %s "
                    "    AND tu.participate = %s "
                    "    AND tv.vector IS NOT NULL", (self.user_id, THREAD_JOINED))

                for d in q2:
                    if not d.user_id == self.user_id:
                        data['thread-users-ids'].append(d.user_id)
                    data['thread-embedding'].append(d.vector[:self.embedding_size])
                if data['thread-embedding']:
                    data['thread-embedding'] = np.array(data['thread-embedding'])
                else:
                    data['thread-embedding'] = np.zeros((1, self.embedding_size), )

                owner_count = Counter(data['thread-users-ids']).most_common()

                # Normalize averaged embedding
                vec_thread = np.sum(data['thread-embedding'], axis=0)
                norm = np.linalg.norm(vec_thread)
                if norm > 0:
                    vec_thread /= norm
            else:
                vec_paper = []
                auth_count = []
                jour_count = []
                vec_thread = []
                owner_count = []

            # Populate fingerprint
            self.embedding = pad_or_trim_vector(vec_paper)
            self.authors_ids = pad_or_trim_vector([i[0] for i in auth_count])
            self.authors_counts = pad_or_trim_vector([i[1] for i in auth_count])
            self.journals_ids = pad_or_trim_vector([i[0] for i in jour_count])
            self.journals_counts = pad_or_trim_vector([i[1] for i in jour_count])
            self.thread_embedding = pad_or_trim_vector(vec_thread)
            self.threads_users_ids = pad_or_trim_vector([i[0] for i in owner_count])
            self.threads_users_counts = pad_or_trim_vector([i[1] for i in owner_count])

            self.save()
            self.set_state('IDL')


