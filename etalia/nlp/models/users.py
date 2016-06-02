# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
import collections

from django.db import models
from django.utils import timezone
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.users.models import UserLibPaper
from etalia.core.models import TimeStampedModel
from etalia.library.models import AuthorPaper
from etalia.core.utils import pad_vector


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

    data = {'ids': [],
            'journal-ids': [],
            'authors-ids': [],
            'date': [],
            'embedding': []
            }

    class Meta:
        unique_together = (('user', 'name'), )

    def __str__(self):
        return '{email}'.format(email=self.user.email)

    @property
    def embedding_size(self):
        return self.model.size

    def update(self):

        data = self.data

        delta = self.user.settings.stream_roll_back_deltatime  # in months
        self.added_after = timezone.now() - timezone.timedelta(days=delta*30.5)

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
                "       AND ulp.date_created >= %s "
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
        auth_count = collections.Counter(
            [auth for auths in data['authors-ids'][start_ind:]
             for auth in auths])
        jour_count = collections.Counter(data['journal-ids'][start_ind:])

        # Normalize averaged embedding
        vec = np.sum(data['embedding'][start_ind:, :], axis=0)
        norm = np.linalg.norm(vec)
        if norm > 0:
            vec /= norm

        # Populate fingerprint
        self.embedding = pad_vector(vec)
        self.authors_ids = pad_vector(list(auth_count.keys()))
        self.authors_counts = pad_vector(list(auth_count.values()))
        self.journals_ids = pad_vector(list(jour_count.keys()))
        self.journals_counts = pad_vector(list(jour_count.values()))

        self.save()


