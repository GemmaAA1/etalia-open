import numpy as np
from scipy.spatial import distance

from django.utils import timezone

from nlp.models import PaperVectors, JournalVectors
from users.models import UserLibPaper


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
        # get_data
        self.seed_data = self.get_data(self.seed_pks)
        self.target_data = self.get_data(self.target_pks)
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
        mat = np.zeros((len(data), self.model.size), dtype=np.float)
        # populate
        for i, entry in enumerate(data):
            mat[i] = np.array(entry['vector'][:self.model.size])
        return mat

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
            created_date = UserLibPaper.objects\
                .filter(
                    userlib=self.user.lib,
                    paper_id__in=self.seed_pks)\
                .values_list('paper__pk',
                             'date_created')

            # Convert date into lapse of time from now in days
            created_date_d = {k: self.convert_date(v)
                                 for k, v in created_date}

            self._created_date_dict = created_date_d
        return self._created_date_dict

    def build_created_date_vec(self, data):
        """Return an array of weights corresponding to created_date ordered
        by paper_pk in data"""
        date_vec = np.zeros(data.count(), dtype=np.float)
        for i, entry in enumerate(data):
            date_vec[i] = self.logist_weight(
                self.created_date_dict[entry['paper_pk']])
        return date_vec

    @property
    def journal_dict(self):
        """Return a dictionary of journal_pk: vector of all journal in seed_data
        """
        if not self._journal_dict:
            # get journal of seed papers
            jpk = [sd['paper__journal__pk'] for sd in self.seed_data]
            jpk = list(set(jpk))
            journal_data = JournalVectors.objects\
                .filter(journal__pk__in=jpk, model=self.model)\
                .values_list('journal__pk', 'vector')
            journal_dict = dict(journal_data)

            self._journal_dict = journal_dict

        return self._journal_dict

    def build_journal_mat(self, data):
        """Return 2D array matching elements of data"""
        journal_mat = np.zeros((data.count(), self.model.size), dtype=np.float)
        for i, entry in enumerate(data):
            journal_mat[i] = np.array(
                self.journal_dict[entry['paper__journal__pk']][:self.model.size])
        return journal_mat

    def weight_with_journal(self, data, mat):
        """Return the weights 2D array of journal and mat"""
        journal_mat = self.build_journal_mat(data)
        return (1. - self.journal_ratio) * mat + self.journal_ratio*journal_mat

    def _run(self):
        """Returns a subset of target primary key and corresponding score"""
        raise NotImplementedError


class SimpleAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        scores = 1. - np.average(distance.cdist(seed_mat, targ_mat, 'cosine'),
                                 axis=0)
        return self.target_pks, scores


class ThresholdAverage(Scoring):

    def __init__(self, **kwargs):
        super(ThresholdAverage, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', 0.6)

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        dis = 1. - distance.cdist(seed_mat, targ_mat, 'cosine')
        dis = np.where(dis > self.threshold, 1, 0)
        scores = np.sum(dis, axis=0)

        return self.target_pks, scores


class WeightedJournalAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        scores = 1. - np.average(distance.cdist(seed_mat, targ_mat, 'cosine'),
                                 axis=0)
        return self.target_pks, scores


class WeightedJournalCreatedDateAverage(Scoring):

    def _run(self):
        seed_mat = self.build_mat(self.seed_data)
        targ_mat = self.build_mat(self.target_data)
        # weight with journal
        seed_mat = self.weight_with_journal(self.seed_data, seed_mat)
        date_vec = self.build_created_date_vec(self.seed_data)

        dis = 1.0 - distance.cdist(seed_mat, targ_mat, 'cosine')
        scores = np.average(dis, weights=date_vec, axis=0)

        return self.target_pks, scores


