# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
from gensim import matutils
from progressbar import ProgressBar, Percentage, Bar, ETA

from django.db import models
from django.db.models import Q
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.core.models import TimeStampedModel
from etalia.core.utils import pad_or_trim_vector, pad_neighbors

from etalia.threads.models import Thread
from etalia.threads.constant import THREAD_TIME_LAPSE_CHOICES
from etalia.library.models import PaperUser
from etalia.library.constants import PAPER_ADDED, PAPER_PINNED
from etalia.users.models import UserSettings

from .users import UserFingerprint


class ThreadVectors(TimeStampedModel):
    """ Table Thread - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use setter/getter functions to set and get vector
    """

    thread = models.ForeignKey(Thread,
                               related_name='vectors')

    model = models.ForeignKey('nlp.Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('thread', 'model')

    def __str__(self):
        return '{thread_pk}/{name}'.format(thread_pk=self.thread.pk,
                                           name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_or_trim_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class ThreadNeighbors(TimeStampedModel):
    """Table for Thread neighbors"""
    # thread
    thread = models.ForeignKey(Thread)

    # Time lapse for neighbors match
    time_lapse = models.IntegerField(default=-1,
                                     choices=THREAD_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    # MostSimilarThread
    te = models.ForeignKey('nlp.ThreadEngine')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{te}/{time_lapse}'.format(
            te=self.te.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.te.model.size]

    class Meta:
        unique_together = ('time_lapse', 'thread', 'te')


class ModelThreadMixin(object):
    """Mixin to Model to add Thread support"""
    THREAD_TOKENIZED_FIELDS = ['title', 'content']

    def infer_thread(self, thread_pk, **kwargs):
        """
        Infer vector for model and thread instance
        """
        vector = self.infer_object(Thread, thread_pk,
                                   self.THREAD_TOKENIZED_FIELDS, **kwargs)

        pv, _ = ThreadVectors.objects.get_or_create(model=self,
                                                    thread_id=thread_pk)

        # store
        pv.set_vector(vector.tolist())

        return thread_pk

    def infer_threads(self, thread_pks, **kwargs):
        """Infer vector for model and all thread in thread_pks list
        """
        # Check inputs
        if isinstance(thread_pks, models.QuerySet):
            thread_pks = list(thread_pks)
        if not isinstance(thread_pks, list):
            raise TypeError(
                'thread_pks must be list or QuerySet, found {0} instead'.format(
                    type(thread_pks)))

        # setup progressbar
        nb_pbar_updates = 100
        nb_threads = len(thread_pks)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=nb_threads, redirect_stderr=True).start()

        for count, thread_pk in enumerate(thread_pks):
            self.infer_thread(thread_pk, **kwargs)
            if not (nb_threads // nb_pbar_updates) == 0:
                if not count % (nb_threads // nb_pbar_updates):
                    pbar.update(count)
        # close progress bar
        pbar.finish()

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented


class ThreadEngineScoringMixin(object):

    SCORE_N_THREADS = 250
    SCORE_USER_CAP_COUNT = 10

    model = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    embedding_size = None

    score_user_boost = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_paper_boost = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_thread_embedding_paper_weight = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    # Matching ThreadEngine class
    data = {'ids': [],
            'paper-ids': [],
            'users-ids': [],
            'date': [],
            'embedding': np.empty(0)
            }

    def score_threadfeed(self, user_id, name='main'):

        # Gather user fingerprint
        if UserFingerprint.objects.filter(user_id=user_id,
                                          name=name,
                                          model=self.model).exists():
            f = UserFingerprint.objects.get(user_id=user_id,
                                            name=name,
                                            model=self.model)
        else:
            f = UserFingerprint.objects.create(user_id=user_id,
                                               name=name,
                                               model=self.model)
            f.update()

        # Get user settings
        us = UserSettings.objects.get(user_id=user_id)

        # Convert user data
        seed = np.array(f.thread_embedding[:self.embedding_size]) \
               + self.score_thread_embedding_paper_weight \
               * np.array(f.embedding[:self.embedding_size])
        ub = self.convert_to_user_boost(f.threads_users_counts)
        ubdic = dict([(k, ub[i]) for i, k in enumerate(f.threads_users_ids)])

        # Compute
        uboost = np.zeros((self.data['embedding'].shape[0], ))
        for i, uid in enumerate(self.data['users-ids']):
            if uid in list(ubdic.keys()):
                uboost[i] = ubdic[uid]

        user_pids = PaperUser.objects.filter(
            Q(user_id=user_id) & (Q(store=PAPER_ADDED) | Q(watch=PAPER_PINNED))) \
            .values_list('paper_id', flat=True)
        pboost = np.zeros((self.data['embedding'].shape[0], ))
        for i, pid in enumerate(self.data['paper-ids']):
            if pid in user_pids:
                pboost[i] = self.score_paper_boost

        # Compute
        score = np.dot(self.data['embedding'], seed.T) + \
                uboost + \
                pboost

        results = self.order_n(self.data['ids'],
                               score,
                               self.data['date'],
                               self.SCORE_N_THREADS)

        return results

    def order_n(self, ids, vals, date, n):
        """"Return top scores

        Return the top N vals and order in decreasing val order
        [{'id': id, 'score': val, 'date': date}, ...]
        """
        # sort scores
        top_idx = matutils.argsort(vals, topn=n, reverse=True)

        return [{'id': ids[i],
                 'score': float(vals[i]),
                 'date': date[i],
                 } for i in top_idx]

    def convert_to_user_boost(self, count):
        return self.convert_to_boost(count,
                                     self.SCORE_USER_CAP_COUNT,
                                     self.score_user_boost)

    def convert_to_boost(self, count, cap, max_boost):
        boost = []
        for c in count:
            if c >= cap:
                boost.append(max_boost)
            else:
                boost.append(c / cap * max_boost)
        return boost