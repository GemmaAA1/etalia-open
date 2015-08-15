from .base import UserFeedTestCase
from ..utils import Scoring, SimpleAverage, ThresholdAverage, \
    WeightedJournalAverage, WeightedJournalCreatedDateAverage
from django.utils import timezone

class ScoringTest(UserFeedTestCase):

    def setUp(self):
        super(ScoringTest, self).setUp()

    def test_scoring_can_be_instantiated(self):
        Scoring(model=self.model, user=self.user)

    def test_scoring_canNOT_be_run_if_not_ready(self):
        scoring = Scoring(model=self.model, user=self.user)
        with self.assertRaises(AssertionError):
            scoring.score()

    def test_get_data_from_pks(self):
        pks = self.papers.values('pk')
        scoring = Scoring(model=self.model, user=self.user)
        data = scoring.get_data(pks)
        self.assertIsNotNone(data)

    def test_get_data_from_pks_in_list(self):
        pks = self.papers.values_list('pk', flat=True)
        scoring = Scoring(model=self.model, user=self.user)
        data = scoring.get_data(pks)
        self.assertIsNotNone(data)

    def test_build_mat_from_data(self):
        pks = self.papers.values('pk')
        scoring = Scoring(model=self.model, user=self.user)
        data = scoring.get_data(pks)
        mat = scoring.build_mat(data)
        self.assertEqual(mat.shape, (pks.count(), self.model.size))

    def test_prepare(self):
        pks = self.papers.values('pk')
        scoring = Scoring(model=self.model, user=self.user)
        scoring.prepare(pks, pks)

    def test_convert_date(self):
        days0 = 10
        date = timezone.datetime.now().date() - timezone.timedelta(days=days0)
        scoring = Scoring(model=self.model, user=self.user)
        days = scoring.convert_date(date)
        self.assertEqual(-days0, days)

    def test_logist_weight_is_within_baseline_and_one(self):
        inputs = [-100, -10, -1, -1e6, 1e6]
        scoring = Scoring(model=self.model, user=self.user)
        vals = [scoring.logist_weight(a) for a in inputs]
        self.assertTrue(all([b >= 0 for b in vals]))
        self.assertTrue(all([b <= 1 for b in vals]))
        self.assertTrue(vals[0] < vals[1] < vals[2])

    def test_created_date_dict(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = Scoring(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        self.assertIsNone(scoring._created_date_dict)
        ddict = scoring.created_date_dict
        self.assertIsNotNone(scoring._created_date_dict)
        self.assertEqual(ddict[self.paper4.pk], -10)

    def test_create_date_vec(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = Scoring(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        vec = scoring.build_created_date_vec

    def test_journal_dict(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = Scoring(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        self.assertIsNone(scoring._journal_dict)
        jdict = scoring.journal_dict
        self.assertIsNotNone(scoring._journal_dict)
        self.assertEqual(jdict[self.journal1.pk],
                         self.journal1.vectors.get(model=self.model)
                         .vector)

    def test_build_journal_mat(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = Scoring(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        j_mat = scoring.build_journal_mat(scoring.seed_data)
        self.assertEqual(j_mat.shape, (len(seed_pks), self.model.size))


class SubScoringTest(UserFeedTestCase):

    def test_simple_average_run(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = SimpleAverage(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        scores, pks = scoring.score()
        self.assertEqual(len(pks), len(targ_pks))
        self.assertEqual(len(scores), len(targ_pks))

    def test_threshold_average_run(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = ThresholdAverage(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        scores, pks = scoring.score()
        self.assertEqual(len(pks), len(targ_pks))
        self.assertEqual(len(scores), len(targ_pks))

    def test_journal_weighted_average_run(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = WeightedJournalAverage(model=self.model, user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        scores, pks = scoring.score()
        self.assertEqual(len(pks), len(targ_pks))
        self.assertEqual(len(scores), len(targ_pks))

    def test_date_and_journal_weighted_average_run(self):
        seed_pks = [p.id for p in self.user.lib.papers.all()]
        targ_pks = [p.id for p in [self.paper, self.paper2]]
        scoring = WeightedJournalCreatedDateAverage(model=self.model,
                                                    user=self.user)
        scoring.prepare(seed_pks, targ_pks)
        scores, pks = scoring.score()
        self.assertEqual(len(pks), len(targ_pks))
        self.assertEqual(len(scores), len(targ_pks))