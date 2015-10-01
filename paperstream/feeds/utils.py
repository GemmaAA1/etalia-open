# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
from scipy.spatial import distance

from django.utils import timezone

from paperstream.nlp.models import PaperVectors, JournalVectors
from paperstream.library.models import Paper

class Scoring(object):
    """Scoring abstract class"""

    def __init__(self, model, user, **kwargs):
        self.model = model
        self.user = user
        self.journal_ratio = kwargs.get('journal_ratio', 0.2)

        self.seed_pks = []
        self.target_pks = []
        self.seed_data = []
        self.target_data = []
        # contains all [journal_pk]:vectors for journals from papers in seed_pks
        self._journal_dict = {}
        self._created_date_dict = {}
        self.is_ready = False

    def prepare(self, seed, target):
        """Score target papers related to seed papers

        Args:
            seed (QuerySet): queryset of Paper primary keys
            target (QuerySet): queryset of Paper primary keys

        """
        # restore to none
        self._journal_dict = None
        self._created_date_dict = None
        # populate
        self.seed_pks = seed
        self.target_pks = target
        # get_data at once
        pks = list(self.seed_pks) + list(self.target_pks)
        data = self.get_data(pks)
        # re-dispatch data
        self.seed_data = [d for d in data if d['paper__pk'] in self.seed_pks]
        self.target_data = [d for d in data if d['paper__pk'] in self.target_pks]
        self.is_ready = True

    def score(self):
        """Score target papers related to seed papers

        Args:
            seed (QuerySet): queryset of Paper primary keys
            target (QuerySet): queryset of Paper primary keys

        Returns:
            (array): array of Paper primary keys (subset of target)
            (array): corresponding scores
        """
        assert self.is_ready
        return self._run()

    def get_data(self, pks):
        data = PaperVectors.objects\
            .filter(
                paper__pk__in=pks,
                model=self.model)\
            .values('vector',
                    'paper__pk',
                    'paper__date_ep',
                    'paper__date_pp',
                    'paper__journal__pk')
        return data

    def build_mat(self, data):
        """Return 2D array of seed paper vectors (papers x vector)"""
        # init seed-mat
        vectors = [d['vector'][:self.model.size] for d in data]
        return np.array(vectors)

    @staticmethod
    def convert_date(date):
        """Convert date into lapse of time from now in days"""
        return (date - timezone.datetime.today().date()).days

    @staticmethod
    def logist_weight(day_lapse, baseline=0.5, k=0.1, delay=-60.0):
        """Logistic function with constant baseline term

        Args:
            day_lapse (int, float): time in days (in the past form now i.e. <0)
            baseline (float): baseline
            k (float): steepness
            delay (float): logistic function kicks off at that time from now

        Returns:
            (float): weighted date by logistic function
        """
        return baseline + (1-baseline) / (1 + np.exp(- (day_lapse - delay) * k))

    @property
    def created_date_dict(self):
        """Return dictionary of paper primary key and created date

        Return dictionary of type {paper_pk: created_date } where created_date
        is in 'time lapse from now in days'"""
        if not self._created_date_dict:
            # Get date data
            created_date = Paper.objects\
                .filter(userlib_paper__userlib=self.user.lib,
                        userlib_paper__paper__in=self.seed_pks)\
                .values_list('pk',
                             'userlib_paper__date_created')

            # Convert date into lapse of time from now in days
            created_date_d = {k: self.convert_date(v) for k, v in created_date}

            self._created_date_dict = created_date_d
        return self._created_date_dict

    def build_created_date_vec(self, data):
        """Return an array of weights corresponding to created_date ordered
        by paper_pk in data"""
        date_vec = np.zeros(data.count(), dtype=np.float)
        for i, entry in enumerate(data):
            date_vec[i] = self.logist_weight(
                self.created_date_dict[entry['paper__pk']])
        return date_vec

    @property
    def journal_dict(self):
        """Return a dictionary of journal_pk: vector of all journal in seed_data
        """
        if not self._journal_dict:
            # get journal of seed papers
            jpk = [sd['paper__journal__pk'] for sd in self.seed_data
                   if sd['paper__journal__pk']]
            jpk = list(set(jpk))
            journal_data = JournalVectors.objects\
                .filter(journal__pk__in=jpk, model=self.model)\
                .values_list('journal__pk', 'vector')
            journal_dict = dict(journal_data)

            self._journal_dict = journal_dict

        return self._journal_dict

    def build_journal_mat(self, data):
        """Return 2D array of journal vectors matching order of papers

        If journal is not defined, set to zeros vector and set corresponding
        ratio to 0
        """
        vectors = []
        zeros_vector = np.zeros((self.model.size, ))
        ratio = []
        for d in data:
            if d['paper__journal__pk']:
                vectors.append(self.journal_dict.get(
                    d['paper__journal__pk'], zeros_vector)[:self.model.size])
                ratio.append(self.journal_ratio)
            else:
                vectors.append(zeros_vector)
                ratio.append(0.0)

        return np.array(vectors), np.array(ratio)

    def weight_with_journal(self, data, mat):
        """Return the weights 2D array of journal and mat"""
        journal_mat, ratio = self.build_journal_mat(data)
        return ((1. - ratio) * mat.T).T + (ratio * journal_mat.T).T

    def _run(self):
        """Returns a subset of target primary key and corresponding score"""
        raise NotImplementedError


class SimpleMax(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        scores = np.max(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class SimpleAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        scores = np.average(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class ThresholdAverage(Scoring):

    def __init__(self, **kwargs):
        super(ThresholdAverage, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', 0.6)

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        dis = np.dot(targ_mat, seed_mat.T)
        dis = np.where(dis > self.threshold, 1, 0)
        scores = np.sum(dis, axis=1)

        return self.target_pks, scores


class WeightedJournalAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        scores = np.average(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class WeightedJournalCreatedDateAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        date_vec = self.build_created_date_vec(self.seed_data)

        dis = np.dot(targ_mat, seed_mat.T)
        scores = np.average(dis, weights=date_vec, axis=1)

        return self.target_pks, scores


