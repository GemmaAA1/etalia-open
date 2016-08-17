# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
from gensim import matutils
from progressbar import ProgressBar, Percentage, Bar, ETA

from django.db import models
from django.conf import settings
from django.contrib.postgres.fields import ArrayField

from etalia.core.models import TimeStampedModel
from etalia.core.utils import pad_or_trim_vector, pad_neighbors

from etalia.library.constants import PAPER_TIME_LAPSE_CHOICES
from etalia.library.models import Paper, Journal
from etalia.users.models import UserSettings

from .users import UserFingerprint


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

    paper = models.ForeignKey(Paper,
                              related_name='vectors')

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
        self.vector = pad_or_trim_vector(vector)
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
        self.vector = pad_or_trim_vector(vector)
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
            pe=self.pe.name,
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


class PaperEngineScoringMixin(object):

    SCORE_AUTHOR_CAP_COUNT = 10
    SCORE_JOURNAL_CAP_COUNT = 10
    SCORE_ALTMETRIC_CAP_SCORE = 200

    score_author_boost = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_journal_boost = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_altmetric_boost = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_n_papers = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_stream_threshold = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    score_trend_threshold = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    model = None  # TO BE ADDED AS MODEL FIELD TO BASE CLASS

    embedding_size = None

    # Matching PaperEngine class
    data = {'ids': [],
            'journal-ids': [],
            'authors-ids': [],
            'date': [],
            'embedding': np.empty(0),
            'altmetric': [],
            }

    def get_user_fingerprint(self, user_id, name='main'):

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

        return f

    def score_stream(self, user_id, name='main'):

        # Gather user fingerprint
        f = self.get_user_fingerprint(user_id, name=name)

        if f.embedding:     # fingerprint is defined (library is not empty)
            # Convert user data
            seed = np.array(f.embedding[:self.embedding_size])
            jb = self.convert_to_journal_boost(f.journals_counts)
            ab = self.convert_to_author_boost(f.authors_counts)
            jbdic = dict([(k, jb[i]) for i, k in enumerate(f.journals_ids)])
            abdic = dict([(k, ab[i]) for i, k in enumerate(f.authors_ids)])

            # Get user settings
            us = UserSettings.objects.get(user_id=user_id)

            # Compute
            jboost = np.zeros((self.data['embedding'].shape[0], ))
            for i, jid in enumerate(self.data['journal-ids']):
                if jid in list(jbdic.keys()):
                    jboost[i] = jbdic[jid]

            aboost = np.zeros((self.data['embedding'].shape[0], ))
            aids_set = set(abdic.keys())
            for i, aids in enumerate(self.data['authors-ids']):
                b = 0
                for aid in aids_set.intersection(set(aids)):
                    b += abdic[aid]
                aboost[i] = min([b, self.score_author_boost])

            score = us.stream_vector_weight * np.dot(self.data['embedding'], seed.T) + \
                    us.stream_journal_weight * jboost + \
                    us.stream_author_weight * aboost

            results = self.order_n(self.data['ids'], score, self.data['date'], self.score_n_papers)

            return results
        return []

    def score_trend(self, user_id, name='main'):

        # Gather user fingerprint
        f = self.get_user_fingerprint(user_id, name=name)

        if f.embedding:     # fingerprint is defined (library is not empty)
            # Convert user data
            seed = np.array(f.embedding[:self.embedding_size])

            # Get user settings
            us = UserSettings.objects.get(user_id=user_id)

            # Convert altmetric
            altmetric_boost = np.array(self.convert_to_boost(
                self.data['altmetric'],
                self.SCORE_ALTMETRIC_CAP_SCORE,
                self.score_altmetric_boost))

            score = us.trend_doc_weight * np.dot(self.data['embedding'], seed.T) + \
                    us.trend_altmetric_weight * altmetric_boost

            results = self.order_n(self.data['ids'], score, self.data['date'], self.score_n_papers)

            return results
        return []

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

    def convert_to_journal_boost(self, count):
        return self.convert_to_boost(count,
                                     self.SCORE_JOURNAL_CAP_COUNT,
                                     self.score_journal_boost)

    def convert_to_author_boost(self, count):
        return self.convert_to_boost(count,
                                     self.SCORE_AUTHOR_CAP_COUNT,
                                     self.score_author_boost)

    def convert_to_boost(self, count, cap, max_boost):
        boost = []
        for c in count:
            if c:
                if c >= cap:
                    boost.append(max_boost)
                else:
                    boost.append(c / cap * max_boost)
            else:
                boost.append(0)
        return boost





