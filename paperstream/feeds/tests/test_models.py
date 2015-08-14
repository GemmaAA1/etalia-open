import numpy as np

from django.core.exceptions import ValidationError

from .base import UserFeedTestCase
from ..models import UserFeed, UserFeedVector, UserFeedPaper


class UserFeedBasicTest(UserFeedTestCase):

    def setUp(self):
        super(UserFeedBasicTest, self).setUp()

    def test_userfeed_can_be_instantiated(self):
        UserFeed(name='test', user=self.user)

    def test_userfeed_must_have_user(self):
        uf = UserFeed()
        with self.assertRaises(ValidationError):
            uf.full_clean()

    def test_userfeed_default_status_is_none(self):
        uf = UserFeed(name='test', user=self.user)
        self.assertEqual(uf.state, 'NON')

    def test_userfeed_default_name_is_main(self):
        uf = UserFeed(user=self.user)
        self.assertEqual(uf.name, 'main')

    def test_userfeed_cannot_have_same_name_and_user(self):
        uf = UserFeed(name='test', user=self.user)
        uf.save()
        uf = UserFeed(name='test', user=self.user)
        with self.assertRaises(ValidationError):
            uf.full_clean()


class UserFeedPaperTest(UserFeedTestCase):

    def setUp(self):
        super(UserFeedPaperTest, self).setUp()
        self.userfeed = UserFeed(name='test', user=self.user)
        self.userfeed.save()

    def test_userfeedpaper_can_be_created(self):
        UserFeedPaper.objects.create(paper=self.paper, feed=self.userfeed)

    def test_userfeedpaper_defaults_are_0_or_none(self):
        ufp = UserFeedPaper(paper=self.paper, feed=self.userfeed)
        self.assertFalse(ufp.is_score_computed)
        self.assertFalse(ufp.is_disliked)
        self.assertFalse(ufp.is_liked)
        self.assertEqual(ufp.score, 0.0)

    def test_userfeedpaper_canNOT_have_same_paper_and_feed(self):
        UserFeedPaper.objects.create(paper=self.paper, feed=self.userfeed)
        ufp = UserFeedPaper(paper=self.paper, feed=self.userfeed)
        with self.assertRaises(ValidationError):
            ufp.full_clean()

    def test_userfeedpaper_are_ordered_by_score(self):
        ufp1 = UserFeedPaper.objects.create(paper=self.paper,
                                           feed=self.userfeed,
                                           score=1.0)
        ufp2 = UserFeedPaper.objects.create(paper=self.paper2,
                                           feed=self.userfeed,
                                           score=2.0)
        ufp3 = UserFeedPaper.objects.create(paper=self.paper3,
                                           feed=self.userfeed,
                                           score=3.0)
        ufps = UserFeedPaper.objects.all()
        self.assertEqual(list(ufps), [ufp3, ufp2, ufp1])


class UserFeedVectorTestCase(UserFeedTestCase):

    def setUp(self):
        super(UserFeedVectorTestCase, self).setUp()
        self.userfeed = UserFeed(name='test', user=self.user)
        self.userfeed.save()

    def test_userfeedvector_can_be_created(self):
        UserFeedVector.objects.create(model=self.model,
                                      feed=self.userfeed)

    def test_userfeedvector_is_unique(self):
        UserFeedVector.objects.create(model=self.model,
                                      feed=self.userfeed)
        ufv = UserFeedVector(model=self.model, feed=self.userfeed)
        with self.assertRaises(ValidationError):
            ufv.full_clean()

    def test_userfeedvector_can_set_vector(self):
        ufv = UserFeedVector.objects.create(model=self.model,
                                            feed=self.userfeed)
        vector = np.random.randn(self.model.size)
        ufv.set_vector(vector)
        self.assertIsNotNone(ufv.vector)

    def test_userfeedvector_can_set_and_get_back_vector(self):
        ufv = UserFeedVector.objects.create(model=self.model,
                                            feed=self.userfeed)
        vector = np.random.randn(self.model.size)
        ufv.set_vector(vector)
        vector2 = ufv.get_vector()
        self.assertEqual(list(vector), vector2)

#
# class UserFeedTest(UserFeedTestCase):
#
#     def setUp(self):
#         super(UserFeedTest, self).setUp()
#         self.userfeed = UserFeed(name='test', user=self.user)
#         self.userfeed.save()
#
#     def test_userfeed_can_add_seed_paper(self):
#         self.userfeed.add_seed_papers(self.papers)
#         self.assertTrue(self.userfeed..papers_seed.count(), self.papers.count())
#
#     def test_userfeed_can_update_seedvector(self):
#         self.userfeed.add_seed_papers(self.papers)
#         self.userfeed.update_userfeed_vector()


    # def test_userfeed_can_be_created(self):
    #     UserFeed.objects.create(name='test', user=self.user)
