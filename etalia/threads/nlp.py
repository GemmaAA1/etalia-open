# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from progressbar import ProgressBar, Percentage, Bar, ETA

from django.db import models
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.core.models import TimeStampedModel
from etalia.core.utils import pad_vector, pad_neighbors

from .constant import THREAD_TIME_LAPSE_CHOICES
from .models import Thread


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
        self.vector = pad_vector(vector)
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
    ms = models.ForeignKey('nlp.MostSimilarThread')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{ms}/{time_lapse}'.format(
            ms=self.ms.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.ms.model.size]

    class Meta:
        unique_together = ('time_lapse', 'thread', 'ms')


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

    def tasks(self, task, **kwargs):
        if task == 'infer_thread':
            thread_pk = kwargs.pop('thread_pk')
            return self.infer_thread(thread_pk, **kwargs)
        if task == 'infer_threads':
            thread_pks = kwargs.pop('thread_pks')
            return self.infer_threads(thread_pks, **kwargs)
        return super(ModelThreadMixin, self).tasks(task, **kwargs)

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented