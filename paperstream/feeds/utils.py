# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
import collections
from scipy.spatial import distance

from django.utils import timezone

from paperstream.nlp.models import PaperVectors, JournalVectors
from paperstream.library.models import Paper


class StreamScoring(object):
    """Scoring abstract class"""

    def __init__(self, model, user, **kwargs):
        self.model = model
        self.user = user
        self.journal_ratio = kwargs.get('journal_ratio', 0.5)
        self.min_auth_cut = kwargs.get('min_auth_cut', 3)
        self.max_auth_cut = kwargs.get('max_auth_cut', 10)
        self.min_auth_p = kwargs.get('min_auth_p', 0.5)
        self.min_jour_cut = kwargs.get('min_jour_cut', 1)
        self.max_jour_cut = kwargs.get('max_jour_cut', 10)
        self.min_jour_p = kwargs.get('min_jour_p', 0.3)
        self.seed_pks = []
        self.seed_auth_pk = []
        self.seed_jour_pk = []
        self.target_pks = []
        self.seed_data = []
        self.target_data = []
        self.seed_auth_data = []
        self.target_auth_data = []
        self.profile = []
        self.vec_w = 1.
        self.auth_w = 0.
        self.jour_w = 0.
        # contains all [journal_pk]:vectors for journals from matches in seed_pks
        self._journal_dict = {}
        self._created_date_dict = {}
        self.is_ready = False

    def prepare(self, seed, target):
        """Score target matches related to seed matches

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
        pks = self.seed_pks + self.target_pks
        data = self.get_data(pks)
        # re-dispatch data and sort
        self.seed_data = [d for d in data if d['paper__pk'] in self.seed_pks]
        self.seed_data.sort(key=lambda d: self.seed_pks.index(d['paper__pk']))
        self.target_data = [d for d in data if d['paper__pk'] in self.target_pks]
        self.target_data.sort(key=lambda d: self.target_pks.index(d['paper__pk']))

        # get authors data at once
        auth_data = self.get_auth_data(pks)
        # re-dispatch data and sort
        self.seed_auth_data = [auth_data[pk] for pk in self.seed_pks]
        self.target_auth_data = [auth_data[pk] for pk in self.target_pks]

        # build seed_auth_pk
        self.build_seed_auth_pk()

        # build seed_jour_pk
        self.build_seed_jour_pk()

        # adjust vec, auth, jour weights
        self.update_weights()

        # build profile
        self.build_profile()



        self.is_ready = True

    def score(self):
        """Score target matches related to seed matches

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

    def get_auth_data(self, pks):
        """Get author data and group them by paper pk"""
        auth_data = Paper.objects\
            .filter(pk__in=pks)\
            .values('pk', 'authors')

        auth_dic = {}
        for x in auth_data:
            if auth_dic.get(x['pk']):
                if x.get('authors'):
                    auth_dic[x['pk']] += [x['authors']]
            else:
                if x.get('authors'):
                    auth_dic[x['pk']] = [x['authors']]

        return auth_dic

    def update_weights(self):
        """Adjust weights based on features length"""
        pass

    def build_mat(self, data):
        """Return 2D array of seed paper vectors (matches x vector)"""
        # init seed-mat
        vectors = [d['vector'][:self.model.size] for d in data]
        return np.array(vectors)

    @staticmethod
    def convert_date(date):
        """Convert date into lapse of time from now in days"""
        return (date - timezone.datetime.today().date()).days

    @staticmethod
    def logist_weight(day_lapse, baseline=0.3, k=0.1, delay=-356.0):
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

    @staticmethod
    def logist_weight_time_score(day_lapse, score, b0=0.5, k0=0.1, d0=-60.0, b1=0, k1=20,
                      d1=0.3):
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
        return (b0 + (1-b0) / (1 + np.exp(- (day_lapse - d0) * k0))) * \
               (b1 + (1-b1) / (1 + np.exp(- (score - d1) * k1)))

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
        date_vec = np.zeros(len(data), dtype=np.float)
        for i, entry in enumerate(data):
            date_vec[i] = self.logist_weight(
                self.created_date_dict[entry['paper__pk']])
        return date_vec

    @property
    def journal_dict(self):
        """Return a dictionary of journal_pk: vector of all journal in seed_data
        """
        if not self._journal_dict:
            # get journal of seed matches
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
        """Return 2D array of journal vectors matching order of matches

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

    def build_p_authors(self):
        """Build dico of probability that author in target paper is interested
        for user"""
        # get authors data
        pks = list(self.seed_pks) + list(self.target_pks)
        data = Paper.objects.filter(pk__in=pks).values('pk', 'authors')
        seed_auth = [d['authors'] for d in data if d['pk'] in self.seed_pks]
        target_auth = [d for d in data if d['pk'] in self.target_pks]
        # reformat target_auth and group by paper
        target_auth_r = {}
        for x in target_auth:
            if target_auth_r.get(x['pk']):
                target_auth_r[x['pk']] += [x['authors']]
            else:
                target_auth_r[x['pk']] = [x['authors']]

        # compute seed_auth probabilities
        # compute occurrences
        seed_dist = collections.Counter(seed_auth)
        # convert to prob
        slope = (1.0 - self.min_auth_p)/(self.max_auth_cut - self.min_auth_cut)
        ord = self.min_auth_p - slope * self.min_auth_cut
        seed_p = {}
        for k, v in seed_dist.items():
            if v >= self.max_auth_cut:
                seed_p[k] = 1.0
            elif v <= self.min_auth_cut:
                seed_p[k] = self.min_auth_p
            else:
                seed_p[k] = slope * v + ord

        # associate target_auth to seed_auth probabilities
        p = {}
        for pk in self.target_pks:
            p[pk] = 0.
            auth_pks = target_auth_r[pk]
            for auth_pk in auth_pks:
                p[pk] = max([seed_p.get(auth_pk, self.min_auth_p), p[pk]])

        return p

    def build_p_journals(self):
        """Build dico of probability that journal in target paper is interested
        for user"""
        # get seed journal data
        seed_jour = [d['paper__journal__pk'] for d in self.seed_data]
        # compute occurrences
        seed_dist = collections.Counter(seed_jour)
        # convert to prob
        slope = (1.0 - self.min_jour_p)/(self.max_jour_cut - self.min_jour_cut)
        ord = self.min_jour_p - slope * self.min_jour_cut
        seed_p = {}
        for k, v in seed_dist.items():
            if v >= self.max_auth_cut:
                seed_p[k] = 1.0
            elif v <= self.min_auth_cut:
                seed_p[k] = self.min_auth_p
            else:
                seed_p[k] = slope * v + ord

        # associate target_auth to seed_auth probabilities
        p = {}
        for pk in self.target_pks:
            p[pk] = seed_p.get(pk, self.min_jour_p)

        return p

    def build_profile(self, time_weight=True):

        # build time array
        if time_weight:
            date_vec = self.build_created_date_vec(self.seed_data)
        else:
            date_vec = np.ones((len(self.seed_data, )))

        # build seed mat
        seed_vec_mat = self.build_mat(self.seed_data)

        # build seed author mat
        seed_auth_mat = self.build_auth_mat(self.seed_auth_data)

        # build journal mat
        seed_jour_mat = self.build_jour_mat(self.seed_data)

        # concatenate
        seed_mat = np.hstack((self.vec_w * seed_vec_mat,
                              self.auth_w * seed_auth_mat,
                              self.jour_w * seed_jour_mat))

        # weight average
        self.profile = np.average(seed_mat, weights=date_vec, axis=0)

        # normalize
        self.profile /= np.linalg.norm(self.profile)

        # return profile
        return self.profile

    def build_auth_mat(self, auth_data):
        # Build author mat
        # build mat
        auth_mat = np.zeros((len(auth_data), len(self.seed_auth_pk)))

        # build auth_mat
        for i, auths in enumerate(auth_data):
            for auth in auths:
                if auth in self.seed_auth_pk:
                    auth_mat[i, self.seed_auth_pk.index(auth)] = 1.
        return auth_mat

    def build_seed_auth_pk(self):
        """Build array of auth pk based on occurrence of author pk in seed
        if occurrence is above threshold"""
        seed_auth = [pk for l in self.seed_auth_data for pk in l]
        # compute occurrences
        seed_occ = collections.Counter(seed_auth)
        # pop key that are below threshold
        for k,v in seed_occ.items():
            if v > self.min_auth_cut:
                self.seed_auth_pk.append(k)

        return self.seed_auth_pk

    def build_jour_mat(self, data):
        jour_data = [d['paper__journal__pk'] for d in data]

        jour_mat = np.zeros((len(jour_data), len(self.seed_jour_pk)))
        # build jour_mat
        for i, jour in enumerate(jour_data):
            if jour in self.seed_jour_pk:
                jour_mat[i, self.seed_jour_pk.index(jour)] = 1.
        return jour_mat

    def build_seed_jour_pk(self):
        """Build array of auth pk based on occurrence of author pk in seed
        if occurrence is above threshold"""
        seed_jour = [d['paper__journal__pk'] for d in self.seed_data if d['paper__journal__pk']]
        # compute occurrences
        seed_occ = collections.Counter(seed_jour)
        # pop key that are below threshold
        for k,v in seed_occ.items():
            if v < self.min_jour_cut:
                seed_occ.pop(k, None)
        # retrieve author pk list
        self.seed_jour_pk = list(seed_occ.keys())

        return self.seed_jour_pk


    def _run(self):
        """Returns a subset of target primary key and corresponding score"""
        raise NotImplementedError


class SimpleMax(StreamScoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        scores = np.max(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class SimpleAverage(StreamScoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        scores = np.average(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class ThresholdAverage(StreamScoring):

    def __init__(self, **kwargs):
        super(ThresholdAverage, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', 0.25)

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        dis = np.dot(targ_mat, seed_mat.T)
        dis = np.where(dis > self.threshold, dis, 0)
        scores = np.sum(dis, axis=1)

        return self.target_pks, scores


class WeightedJournalAverage(StreamScoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        targ_mat = self.weight_with_journal(self.target_data, targ_mat)
        scores = np.max(np.dot(targ_mat, seed_mat.T), axis=1)
        return self.target_pks, scores


class WeightedJournalCreatedDateAverage(StreamScoring):

    def __init__(self, **kwargs):
        super(WeightedJournalCreatedDateAverage, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', 0.25)

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        targ_mat = self.weight_with_journal(self.target_data, targ_mat)
        date_vec = self.build_created_date_vec(self.seed_data)

        dis = np.dot(targ_mat, seed_mat.T)
        dis = np.where(dis > self.threshold, dis, 0)
        scores = np.average(dis, weights=date_vec, axis=1)

        return self.target_pks, scores


class ContentBasedProfile(StreamScoring):

    def __init__(self, **kwargs):
        super(ContentBasedProfile, self).__init__(**kwargs)
        user = kwargs.get('user')

    def _run(self):

        # build target mat
        target_vec_mat = self.build_mat(self.target_data)

        # build target author mat
        target_auth_mat = self.build_auth_mat(self.target_auth_data)

        # build target journal mat
        target_jour_mat = self.build_jour_mat(self.target_data)

        # concatenate
        target_mat = np.hstack((self.vec_w * target_vec_mat,
                                self.auth_w * target_auth_mat,
                                self.jour_w * target_jour_mat))

        # normalize
        target_mat /= np.linalg.norm(target_mat, axis=1)[:, None]

        # dot product
        dis = np.dot(target_mat, self.profile.T)

        return self.target_pks, dis



class OccurrenceCount(StreamScoring):

    def __init__(self, **kwargs):
        super(OccurrenceCount, self).__init__(**kwargs)
        user = kwargs.get('user')
        if user:
            # number of closest neighbors to keep per seed paper
            self.cutoff = 5 + 3 * user.settings.stream_narrowness

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        targ_mat = self.weight_with_journal(self.target_data, targ_mat)
        # build date vector
        # date_vec = self.build_created_date_vec(self.seed_data)
        date_vec = np.ones((len(self.seed_data), ))

        dis = 1.0 - np.dot(targ_mat, seed_mat.T)
        ind = np.argpartition(dis, self.cutoff, axis=0)[:self.cutoff, :][:]

        ind_unique = np.unique(ind[:])

        # count occurrences
        occ = []
        # replicate date_vec for computation ease
        date_mat = np.tile(date_vec, (ind.shape[0], 1))
        for i, idx in enumerate(ind_unique):
            occ.append((idx, np.sum(date_mat[ind == idx])))
        # normalize by number of paper in user lib
        occ = list(map(lambda x: (x[0], x[1]/(seed_mat.shape[0]*self.cutoff)), occ))
        occ_sorted = sorted(occ, key=lambda x: x[1], reverse=True)

        scores = [x[1] for x in occ_sorted]
        scores_pks = [self.target_pks[x[0]] for x in occ_sorted]

        return scores_pks, scores




class TrendScoring(object):
    pass



