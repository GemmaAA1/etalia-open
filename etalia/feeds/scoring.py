# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
import collections
from config.celery import celery_app as app
from scipy.spatial import distance
from gensim import matutils

from django.utils import timezone
from django.conf import settings
from django.core.cache import caches

from etalia.nlp.models import PaperVectors, Model
from etalia.library.models import Paper

from etalia.altmetric.models import AltmetricModel


class Scoring(object):
    """Scoring abstract class"""

    cache = caches['default']
    cache_files = caches['files']

    def __init__(self, stream, **kwargs):
        self.stream = stream
        self.journal_ratio = kwargs.get('journal_ratio', 0.)

        self.is_ready = False
        self.model = Model.objects.get(is_active=True)

        # Init seed data
        self.seed_pks = self.get_seed_pks()
        self.seed_data = self.get_data(self.seed_pks)
        self.seed_auth_data = self.get_auth_data(self.seed_pks)

        # Target data
        self.target_pks = None
        self.target_date = None
        self.target_data = None
        self.target_auth_data = None

        self.journal_ind2jourpk = []
        self.created_date_vec = []

    def get_seed_pks(self):
        """Retrieve seed paper pk"""
        # TODO: optimze SQL query

        seed_pks = PaperVectors.objects \
            .filter(paper__pk__in=self.stream.seeds.values('pk'),
                    model=self.stream.user.settings.stream_model) \
            .values_list('paper__pk', flat=True)

        return list(seed_pks)

    def get_target_all_pks(self):
        """Get all recent paper pks from MostSimilar object"""
        # Get most similar instance
        try:
            ms_task = app.tasks['etalia.nlp.tasks.mostsimilar_{name}'.format(
                name=self.stream.user.settings.stream_model.name)]
        except KeyError:
            raise KeyError
        # Get last pks
        res = ms_task.delay('get_recent_pks',
                            time_lapse=self.stream.user.settings.stream_time_lapse)
        # wait for results
        target_pks, pks_date = res.get()

        # Remove paper already in library
        lib_pks = self.stream.user.lib.papers.all().values_list('pk', flat=True)
        target_pks = [pk for pk in target_pks
                      if pk not in list(lib_pks)][:settings.FEED_MAX_TARGETS]
        pks_date = dict([(pk, pks_date[i]) for i, pk in enumerate(target_pks)
                         if pk not in list(lib_pks)][
                        :settings.FEED_MAX_TARGETS])

        return target_pks, pks_date

    def get_target_neigh_pks_from_seed(self, seed):
        """Get all recent paper pks in the neighborhood of seed using
        MostSimilar"""
        # Get most similar instance
        try:
            ms_task = app.tasks['etalia.nlp.tasks.mostsimilar_{name}'.format(
                name=self.stream.user.settings.stream_model.name)]
        except KeyError:
            raise KeyError

        # Get last pks
        res = ms_task.delay('knn_search',
                            seed=seed,
                            time_lapse=self.stream.user.settings.stream_time_lapse,
                            top_n=settings.FEED_MAX_TARGETS)

        # wait for results
        target_pks = res.get()

        # Remove paper already in library
        lib_pks = self.stream.user.lib.papers.all().values_list('pk', flat=True)
        target_pks = [pk for pk, d in target_pks
                      if pk not in list(lib_pks)][:settings.FEED_MAX_TARGETS]

        return target_pks

    def get_data(self, pks):
        """Get paper metadata and vector and sort by pks"""
        data = list(
            PaperVectors.objects \
                .filter(paper__pk__in=pks, model=self.model) \
                .values('vector', 'paper__pk', 'paper__journal__pk')
        )
        # sort
        data.sort(key=lambda d: pks.index(d['paper__pk']))

        return data

    @staticmethod
    def get_altmetric_data(pks):
        """Get altmetric data"""
        data = list(
            AltmetricModel.objects \
                .filter(paper_id__in=pks) \
                .values('paper__pk', 'score'))
        # sort
        data.sort(key=lambda d: pks.index(d['paper__pk']))

        return data

    def get_auth_data(self, pks):
        """Get author data and group them by paper_pk and sort by pks"""
        papers = Paper.objects\
            .filter(pk__in=pks)\
            .prefetch_related('authors')
        mapping = dict((p.pk, p.authors.all()) for p in papers)
        data = [mapping[x] for x in pks]
        return data

    def build_paper_vec_mat(self, data):
        """Build 2D array of paper vectors (matches x vector)"""
        # init seed-mat
        vectors = [d['vector'][:self.model.size] for d in data]
        return np.array(vectors)

    @staticmethod
    def convert_date(date):
        """Convert date into lapse of time from now in days"""
        return (date - timezone.datetime.today().date()).days

    @staticmethod
    def logist_weight(day_lapse, baseline=0.3, k=0.1, delay=-180):
        """Logistic function with constant baseline term

        Args:
            day_lapse (int, float): time in days (in the past form now i.e. <0)
            baseline (float): baseline
            k (float): steepness
            delay (float): logistic function kicks off at that time from now

        Returns:
            (float): weighted date by logistic function
        """
        return baseline + (1 - baseline) / (
        1 + np.exp(- (day_lapse - delay) * k))

    @staticmethod
    def logist_weight_time_score(day_lapse, score, b0=0.5, k0=0.1, d0=-60.0,
                                 b1=0, k1=20, d1=0.3):
        """2D Logistic (product of 2 1d-logit)

        Args:
            day_lapse (int, float): time in days (in the past form now i.e. <0)
            b0 (float): baseline
            k0 (float): steepness
            d0 (float): logistic function kicks off at that time from now
            score (int, float): score / distance
            b1 (float): baseline
            k1 (float): steepness
            d1 (float): logistic function kicks off at that time from now

        Returns:
            (float): weighted date by logistic function
        """
        return (b0 + (1 - b0) / (1 + np.exp(- (day_lapse - d0) * k0))) * \
               (b1 + (1 - b1) / (1 + np.exp(- (score - d1) * k1)))

    def build_date_vec(self):
        """Build created_date_vec an array of weights related to created_date
        ordered by seed_pks"""
        # init
        date_vec = np.zeros(len(self.seed_pks), dtype=np.float)

        # Get date data
        created_date = self.stream.user.lib.userlib_paper \
            .filter(paper_id__in=self.seed_pks) \
            .values_list('paper_id', 'date_created')

        # Convert date into integer number of days from now
        created_date_int = [(pk, self.convert_date(v)) for pk, v in
                            created_date]

        # logist parameter
        delay = - self.stream.user.settings.stream_roll_back_deltatime * 30
        k = 0.1
        baseline = 0
        created_date_int_d = dict(created_date_int)
        for pk in self.seed_pks:
            w = self.logist_weight(created_date_int_d[pk],
                                   delay=delay,
                                   k=k,
                                   baseline=baseline)
            self.created_date_vec.append(w)

        return self.created_date_vec

    def order_and_trim_results(self, pks, scores):
        """"Sort scores in descending order and trim results
        """
        # number of paper to keep
        nb_papers = int(
            settings.FEED_SIZE_PER_DAY *
            self.stream.user.settings.stream_time_lapse)

        # sort scores
        best = matutils.argsort(scores,
                                topn=nb_papers,
                                reverse=True)

        return dict([(pks[ind], float(scores[ind])) for ind in best])

    def build_doc_profile(self, time_weight=True):

        # build time array
        if time_weight:
            date_vec = self.build_date_vec()
        else:
            date_vec = np.ones((len(self.seed_data, )))

        # build seed mat
        seed_mat = self.build_paper_vec_mat(self.seed_data)

        # weight average
        profile = np.average(seed_mat, weights=date_vec, axis=0)

        # normalize
        norm = np.linalg.norm(profile)
        if norm > 0:
            profile /= norm

        # return profile
        return profile

    def score(self):
        """Score target matches related to seed matches if ready
        """
        raise NotImplemented


class ContentBasedScoring(Scoring):
    """Content-Based Recommendation (CB)

    Match between user profile feature array and target paper feature array.

    Feature array is defined as follow:
     [v1, .... vN, auth1, .... authP, jour1, ... jourQ]

     where:
        v1, ... vN are the feature vector from Doc2Vec
        auth1, ... authP are most cited authors in seeds
        jour1, ... jourQ are most cited journal in seeds

      User profile array is defined as the weighted average of the seed feature
      array

    """

    MAX_SEED_AUTH_PK = 300
    MAX_SEED_JOUR_PK = 100
    MIN_OCC_SEED_AUTH = 2
    MIN_OCC_SEED_JOUR = 1
    DEFAULT_VECTOR_WEIGHT = 1.
    DEFAULT_AUTHOR_WEIGHT = 1.
    DEFAULT_JOURNAL_WEIGHT = 1.
    DEFAULT_TARGET_SEARCH = 'all'

    def __init__(self, **kwargs):
        super(ContentBasedScoring, self).__init__(**kwargs)
        self.vec_w = kwargs.get('vector_weight', self.DEFAULT_VECTOR_WEIGHT)
        self.auth_w = kwargs.get('author_weight', self.DEFAULT_AUTHOR_WEIGHT)
        self.jour_w = kwargs.get('journal_weight', self.DEFAULT_JOURNAL_WEIGHT)

        self.target_search = kwargs.get('target_search',
                                        self.DEFAULT_TARGET_SEARCH)

        # seed index to author_pk
        self.seed_ind2authpk = []
        # seed index to journal_pk
        self.seed_ind2jourpk = []
        # profile array
        self.profile = []

    def build_profile_ind2authpk(self):
        """Build array of author_pk based on occurrence of author_pk in seeds"""
        seed_auth = [pk for l in self.seed_auth_data for pk in l]
        # compute occurrences
        seed_occ = collections.Counter(seed_auth)
        # populate seed_ind2authpk
        self.seed_ind2authpk = []
        for k, v in seed_occ.items():
            if v > self.MIN_OCC_SEED_AUTH:
                self.seed_ind2authpk.append(k)
            if len(self.seed_ind2authpk) >= self.MAX_SEED_AUTH_PK:
                break

        return None

    def build_profile_ind2jourpk(self):
        """Build array of journal_pk based on occurrence of journal_pk in seeds"""
        seed_jour = [d['paper__journal__pk'] for d in self.seed_data if
                     d['paper__journal__pk']]
        # compute occurrences
        seed_occ = collections.Counter(seed_jour)
        self.seed_ind2jourpk = []
        for k, v in seed_occ.items():
            if v > self.MIN_OCC_SEED_JOUR:
                self.seed_ind2jourpk.append(k)
            if len(self.seed_ind2jourpk) >= self.MAX_SEED_JOUR_PK:
                break

        return None

    def build_auth_utility_mat(self, auth_data):
        """Build author utility matrix"""
        auth_mat = np.zeros((len(auth_data), len(self.seed_ind2authpk)))

        # build auth_mat
        for i, auths in enumerate(auth_data):
            for auth in auths:
                if auth in self.seed_ind2authpk:
                    auth_mat[i, self.seed_ind2authpk.index(auth)] = 1.
        return auth_mat

    def build_jour_utility_mat(self, data):
        """Build journal utility matrix"""
        jour_data = [d['paper__journal__pk'] for d in data]

        jour_mat = np.zeros((len(jour_data), len(self.seed_ind2jourpk)))
        # build jour_mat
        for i, jour in enumerate(jour_data):
            if jour in self.seed_ind2jourpk:
                jour_mat[i, self.seed_ind2jourpk.index(jour)] = 1.
        return jour_mat

    def build_profile(self, time_weight=True):

        # filter relevant author
        self.build_profile_ind2authpk()
        # filter relevant journal
        self.build_profile_ind2jourpk()

        # build time array
        if time_weight:
            date_vec = self.build_date_vec()
        else:
            date_vec = np.ones((len(self.seed_data, )))

        # build seed mat
        seed_vec_mat = self.build_paper_vec_mat(self.seed_data)

        # build seed author mat
        seed_auth_mat = self.build_auth_utility_mat(self.seed_auth_data)

        # build journal mat
        seed_jour_mat = self.build_jour_utility_mat(self.seed_data)

        # concatenate these 3 mats
        # seed_mat = np.hstack((self.vec_w * seed_vec_mat,
        #                       self.auth_w * seed_auth_mat,
        #                       self.jour_w * seed_jour_mat))

        # average with time-stamps logistic weights
        seed_vec_av = np.average(seed_vec_mat, weights=date_vec, axis=0)
        seed_auth_av = np.average(seed_auth_mat, weights=date_vec, axis=0)
        seed_jour_av = np.average(seed_jour_mat, weights=date_vec, axis=0)

        # Squashing journal using softmax
        seed_jour_av = np.exp(seed_jour_av) / np.sum(np.exp(seed_jour_av))
        # Squashing authors using softmax
        seed_auth_av = np.exp(seed_auth_av) / np.sum(np.exp(seed_auth_av))

        # concatenate with weight
        self.profile = np.hstack((self.vec_w * seed_vec_av,
                                  self.auth_w * seed_auth_av *
                                  seed_vec_av.shape[0],
                                  self.jour_w * seed_jour_av *
                                  seed_vec_av.shape[0]))

        # # normalize
        # norm = np.linalg.norm(seed_av, axis=0)
        # non_zeros = norm > 0.
        # seed_mat[non_zeros, :] /= norm[non_zeros, None]

        # # weight average
        # self.profile = np.average(seed_mat, weights=date_vec, axis=0)

        # normalize
        norm = np.linalg.norm(self.profile)
        if norm > 0:
            self.profile /= norm

        # return profile
        return self.profile

    def score(self):

        # Build profile
        # NB: should not be be cached (user settings specific)
        self.build_profile(time_weight=True)

        # Get seed
        seed = self.profile[:self.model.size].copy()
        norm = np.linalg.norm(seed)
        if norm > 0:
            seed /= norm

        # Get target data (from cache if available)
        if 'target_pks' not in self.cache:
            if self.target_search == 'neighbor':
                self.target_pks = self.get_target_neigh_pks_from_seed(seed=seed)
            elif self.target_search == 'all':
                self.target_pks, self.target_date = self.get_target_all_pks()
            else:
                raise ValueError('')
            self.cache.add('target_pks', self.target_pks)
            self.cache.add('target_date', self.target_date)
        else:
            self.target_pks = self.cache.get('target_pks')
            self.target_date = self.cache.get('target_date')

        # build target mat (or get from cache)
        user_target_mat_key = '{user_pk}_stream_target_mat'.format(
            user_pk=self.stream.user_id)
        if user_target_mat_key not in self.cache_files:

            # build target data
            if 'target_data' not in self.cache:
                self.target_data = self.get_data(self.target_pks)
                self.cache.add('target_data', self.target_data)
            else:
                self.target_data = self.cache.get('target_data')
            # auth data
            if 'target_auth_data' not in self.cache:
                self.target_auth_data = self.get_auth_data(self.target_pks)
                self.cache.add('target_auth_data', self.target_auth_data)
            else:
                self.target_auth_data = self.cache.get('target_auth_data')

            target_vec_mat = self.build_paper_vec_mat(self.target_data)
            # build target author mat
            target_auth_mat = self.build_auth_utility_mat(self.target_auth_data)
            # build target journal mat
            target_jour_mat = self.build_jour_utility_mat(self.target_data)
            # concatenate
            # target_mat = np.hstack((self.vec_w * target_vec_mat,
            #                         self.auth_w * target_auth_mat,
            #                         self.jour_w * target_jour_mat))
            target_mat = np.hstack((target_vec_mat,
                                    target_auth_mat,
                                    target_jour_mat))

            # normalize
            norm = np.linalg.norm(target_mat, axis=1)
            non_zeros = norm > 0.
            target_mat[non_zeros, :] /= norm[non_zeros, None]
            # cache for user specific target mat
            self.cache_files.add(user_target_mat_key, target_mat)
        else:
            target_mat = self.cache_files.get(user_target_mat_key)

        # Dot product
        dis = np.dot(target_mat, self.profile.T)

        # Order
        res = self.order_and_trim_results(self.target_pks, dis)

        return res, self.target_date


class TrendScoring(Scoring):
    DEFAULT_DOC_WEIGHT = 1.
    DEFAULT_ALTMETRIC_WEIGHT = 1.
    DEFAULT_TARGET_SEARCH = 'all'

    def __init__(self, **kwargs):
        super(TrendScoring, self).__init__(**kwargs)
        self.doc_w = kwargs.get('doc_weight', self.DEFAULT_DOC_WEIGHT)
        self.alt_w = kwargs.get('altmetric_weight',
                                self.DEFAULT_ALTMETRIC_WEIGHT)

        self.target_search = kwargs.get('target_search',
                                        self.DEFAULT_TARGET_SEARCH)

        # profile array
        self.profile = []

    def score(self):

        # Build profile
        self.profile = self.build_doc_profile(time_weight=True)

        # Get seed
        seed = self.profile[:self.model.size].copy()
        norm = np.linalg.norm(seed)
        if norm > 0:
            seed /= norm

        # Get target data (from cache if available)
        # build target mat (or get from cache)
        if 'target_pks' not in self.cache:
            if self.target_search == 'neighbor':
                self.target_pks = self.get_target_neigh_pks_from_seed(seed=seed)
            elif self.target_search == 'all':
                self.target_pks, self.target_date = self.get_target_all_pks()
            else:
                raise ValueError('')
            self.cache.add('target_pks', self.target_pks)
            self.cache.add('target_date', self.target_date)
        else:
            self.target_pks = self.cache.get('target_pks')
            self.target_date = self.cache.get('target_date')

        user_target_mat_key = '{user_pk}_trend_target_mat'.format(
            user_pk=self.stream.user_id)
        if user_target_mat_key not in self.cache_files:
            # Build target data (from cache if available)
            if 'target_data' not in self.cache:
                self.target_data = self.get_data(self.target_pks)
                self.cache.add('target_data', self.target_data)
            else:
                self.target_data = self.cache.get('target_data')
            # auth data
            if 'target_auth_data' not in self.cache:
                self.target_auth_data = self.get_auth_data(self.target_pks)
                self.cache.add('target_auth_data', self.target_auth_data)
            else:
                self.target_auth_data = self.cache.get('target_auth_data')

            # build target mat
            target_mat = self.build_paper_vec_mat(self.target_data)

            # normalize
            norm = np.linalg.norm(target_mat, axis=1)
            non_zeros = norm > 0.
            target_mat[non_zeros, :] /= norm[non_zeros, None]
            # cache for user specific target mat
            self.cache_files.add(user_target_mat_key, target_mat)
        else:
            target_mat = self.cache_files.get(user_target_mat_key)

        # Dot product
        dis = np.dot(target_mat, self.profile.T)

        # Get Altmetric data (from cache if available)
        if 'target_altmetric' not in self.cache:
            target_altmetric = self.get_altmetric_data(self.target_pks)
            # cache
            self.cache.add('target_altmetric', target_altmetric)
        else:
            # get from cache
            target_altmetric = self.cache.get('target_altmetric')
        # Get altmetric score and normalize
        alt_score = np.array([p['score'] for p in target_altmetric])
        # normalize in [0,1]
        alt_score /= alt_score.max()
        # get altmetric pks
        alt_pks = [p['paper__pk'] for p in target_altmetric]

        # Intersect altmetric pks and target pks
        intersect = set.intersection(set(self.target_pks), set(alt_pks))
        # reduced distance based on intersect and get corresponding pks
        pk_idx = np.array([[pk in intersect, pk] for pk in self.target_pks])
        dis_reduced = dis[pk_idx[:, 0] == 1]
        target_pks_reduced = pk_idx[pk_idx[:, 0] == 1, 1]

        # compute weighted score
        score = self.doc_w * dis_reduced + self.alt_w * alt_score

        # Order
        res = self.order_and_trim_results(target_pks_reduced, score)

        return res, self.target_date
