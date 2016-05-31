# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from progressbar import ProgressBar, Percentage, Bar, ETA

from django.db import models
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.core.models import TimeStampedModel
from etalia.core.utils import pad_vector, pad_neighbors

from etalia.library.constants import PAPER_TIME_LAPSE_CHOICES
from etalia.library.models import Paper, Journal


# --------------------------------------
# Library NLP
# --------------------------------------
class PaperVectors(TimeStampedModel):
    """ Table Paper - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use set_vector() to pad and set vector list
    """

    paper = models.ForeignKey(Paper, related_name='vectors')

    model = models.ForeignKey('nlp.Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('paper', 'model')

    def __str__(self):
        return '{paper_pk}/{name}'.format(paper_pk=self.paper.pk,
                                          name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class JournalVectors(TimeStampedModel):
    """ Table Journal - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None.

    Use set_vector() to pad and set vector list
    """

    journal = models.ForeignKey(Journal, related_name='vectors')

    model = models.ForeignKey('nlp.Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('journal', 'model')

    def __str__(self):
        return '{short_title}/{name}' \
            .format(short_title=self.journal.short_title, name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class PaperNeighbors(TimeStampedModel):
    """ Table of matches nearest neighbors"""

    paper = models.ForeignKey(Paper, related_name='neighbors')

    pe = models.ForeignKey('nlp.PaperEngine')

    time_lapse = models.IntegerField(default=-1,
                                     choices=PAPER_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{pe}/{time_lapse}'.format(
            ms=self.pe.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.pe.model.size]

    class Meta:
        unique_together = ('time_lapse', 'paper', 'pe')


class JournalNeighbors(TimeStampedModel):
    """ Table of matches nearest neighbors"""

    journal = models.ForeignKey(Journal, related_name='neighbors')

    pe = models.ForeignKey('nlp.PaperEngine')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{pe}/{time_lapse}'.format(
            pe=self.pe.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.pe.model.size]

    class Meta:
        unique_together = ('journal', 'pe')


class ModelLibraryMixin(object):
    """Mixin to Model to add Thread support"""

    PAPER_TOKENIZED_FIELDS = ['title', 'abstract']

    def infer_paper(self, paper_pk, **kwargs):
        """
        Infer vector for model and thread instance
        """
        vector = self.infer_object(Paper, paper_pk,
                                   self.PAPER_TOKENIZED_FIELDS, **kwargs)

        pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                   paper_id=paper_pk)

        # store
        pv.set_vector(vector.tolist())

        return paper_pk

    def infer_papers(self, paper_pks, **kwargs):
        """Infer vector for model and all thread in thread_pks list
        """
        # Check inputs
        if isinstance(paper_pks, models.QuerySet):
            paper_pks = list(paper_pks)
        if not isinstance(paper_pks, list):
            raise TypeError(
                'thread_pks must be list or QuerySet, found {0} instead'.format(
                    type(paper_pks)))

        # setup progressbar
        nb_pbar_updates = 100
        nb_papers = len(paper_pks)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=nb_papers, redirect_stderr=True).start()

        for count, paper_pk in enumerate(paper_pks):
            self.infer_paper(paper_pk, **kwargs)
            if not (nb_papers // nb_pbar_updates) == 0:
                if not count % (nb_papers // nb_pbar_updates):
                    pbar.update(count)
        # close progress bar
        pbar.finish()

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented