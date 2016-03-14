# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np

from django.core.exceptions import ValidationError
from django.utils import timezone

from .base import UserFeedTestCase
from ..models import Stream, StreamVector, StreamMatches, \
    StreamSeeds
from etalia.library.models import Paper
from etalia.nlp.models import PaperVectors


class UserFeedBasicTest(UserFeedTestCase):

    def test_userfeed_can_be_instantiated(self):
        Stream(name='test', user=self.user)

    def test_userfeed_must_have_user(self):
        uf = Stream()
        with self.assertRaises(ValidationError):
            uf.full_clean()

    def test_userfeed_default_status_is_none(self):
        uf = Stream(name='test', user=self.user)
        self.assertEqual(uf.state, 'NON')

    def test_userfeed_default_name_is_main(self):
        uf = Stream(user=self.user)
        self.assertEqual(uf.name, 'main')

    def test_userfeed_cannot_have_same_name_and_user(self):
        uf = Stream(name='test', user=self.user)
        uf.save()
        uf = Stream(name='test', user=self.user)
        with self.assertRaises(ValidationError):
            uf.full_clean()

    def test_userfeed_can_add_seed_paper(self):
        uf = Stream(name='test', user=self.user)
        uf.save()
        uf.add_papers_seed(self.papers)
        self.assertTrue(uf.seeds.count(), self.papers.count())


class UserFeedPaperTest(UserFeedTestCase):

    def setUp(self):
        super(UserFeedPaperTest, self).setUp()
        self.userfeed = Stream(name='test', user=self.user)
        self.userfeed.save()

    def test_userfeedpaper_can_be_created(self):
        StreamMatches.objects.create(paper=self.paper, feed=self.userfeed)

    def test_userfeedpaper_defaults_are_0_or_none(self):
        ufp = StreamMatches(paper=self.paper, feed=self.userfeed)
        self.assertFalse(ufp.is_score_computed)
        self.assertFalse(ufp.is_disliked)
        self.assertFalse(ufp.is_liked)
        self.assertEqual(ufp.score, 0.0)

    def test_userfeedpaper_canNOT_have_same_paper_and_feed(self):
        StreamMatches.objects.create(paper=self.paper, feed=self.userfeed)
        ufp = StreamMatches(paper=self.paper, feed=self.userfeed)
        with self.assertRaises(ValidationError):
            ufp.full_clean()

    def test_userfeedpaper_are_ordered_by_score(self):
        ufp1 = StreamMatches.objects.create(paper=self.paper,
                                           feed=self.userfeed,
                                           score=1.0)
        ufp2 = StreamMatches.objects.create(paper=self.paper2,
                                           feed=self.userfeed,
                                           score=2.0)
        ufp3 = StreamMatches.objects.create(paper=self.paper3,
                                           feed=self.userfeed,
                                           score=3.0)
        ufps = StreamMatches.objects.all()
        self.assertEqual(list(ufps), [ufp3, ufp2, ufp1])


class UserFeedVectorTestCase(UserFeedTestCase):

    def setUp(self):
        super(UserFeedVectorTestCase, self).setUp()
        self.userfeed = Stream(name='test', user=self.user)
        self.userfeed.save()

    def test_userfeedvector_can_be_created(self):
        StreamVector.objects.create(model=self.model,
                                      feed=self.userfeed)

    def test_userfeedvector_is_unique(self):
        StreamVector.objects.create(model=self.model,
                                      feed=self.userfeed)
        ufv = StreamVector(model=self.model, feed=self.userfeed)
        with self.assertRaises(ValidationError):
            ufv.full_clean()

    def test_userfeedvector_can_set_vector(self):
        ufv = StreamVector.objects.create(model=self.model,
                                            feed=self.userfeed)
        vector = np.random.randn(self.model.size)
        ufv.set_vector(vector)
        self.assertIsNotNone(ufv.vector)

    def test_userfeedvector_can_set_and_get_back_vector(self):
        ufv = StreamVector.objects.create(model=self.model,
                                            feed=self.userfeed)
        vector = np.random.randn(self.model.size)
        ufv.set_vector(vector)
        vector2 = ufv.get_vector()
        self.assertEqual(list(vector), vector2)

    def test_userfeedvector_vector_default_is_zero(self):
        uf = StreamVector.objects.create(model=self.model,
                                           feed=self.userfeed)
        vector = uf.get_vector()
        self.assertEqual(vector, None)

    def test_userfeedvector_can_update_vector(self):
        StreamSeeds.objects.create(feed=self.userfeed,
                                         paper=self.paper)
        StreamSeeds.objects.create(feed=self.userfeed,
                                         paper=self.paper2)
        ufv = StreamVector.objects.create(model=self.model,
                                            feed=self.userfeed)
        ufv.update_vector()


class UserFeedTest(UserFeedTestCase):

    def setUp(self):
        super(UserFeedTest, self).setUp()
        self.userfeed = Stream(name='test', user=self.user)
        self.userfeed.save()

    def test_userfeed_can_add_paper_and_update_seedvector(self):
        self.userfeed.add_papers_seed(self.papers)
        self.assertTrue(self.userfeed.seeds.count() > 0)
        self.assertTrue(self.userfeed.vectors.count() > 0)

    def test_userfeed_can_update(self):
        self.userfeed.add_papers_seed(self.user.lib.papers.all())
        self.user.settings.scoring_method = 1
        self.assertEqual(self.userfeed.matches.count(), 0)
        self.userfeed.update()
        self.assertTrue(self.userfeed.matches.count() > 0)

    def test_userfeed_can_be_created(self):
        ul = Stream.objects.create(user=self.user, name='test2',
                                     papers_seed=self.papers)
        self.assertTrue(ul.name == 'test2')
        self.assertTrue(ul.papers_seed.count() > 0)

    def test_userfeed_can_create_default(self):
        ul = Stream.objects.create_main(user=self.user)
        self.assertEqual(ul.name, 'main')

    def test_userfeed_can_create_and_then_update(self):
        ul = Stream.objects.create(user=self.user)
        ul.add_papers_seed(self.papers)
        ul.update()
        count1 = ul.papers_match.count()
        # a new paper is coming
        new_paper = Paper.objects.create(
            title='Bla bla bla.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=True)
        pv = PaperVectors.objects.create(paper=new_paper, model=self.model)
        vec = np.random.randn(self.model.size)
        pv.set_vector(vec)
        self.model.build_most_similar()
        ul.update()
        count2 = ul.papers_match.count()
        self.assertTrue(count1 + 1 == count2)



