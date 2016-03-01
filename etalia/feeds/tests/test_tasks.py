# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from ..tasks import update_stream, init_main
from ..models import Stream

from .base import UserFeedTestCase

class UpdateTest(UserFeedTestCase):

    def setUp(self):
        super(UpdateTest, self).setUp()
        self.userfeed = Stream.objects.create_main(user=self.user)
        self.userfeed.add_papers_seed(self.papers)

    def test_update_async(self):
        res = update_stream.delay(self.userfeed.pk)
        self.assertTrue(res.successful())


class InitMainFeed(UserFeedTestCase):

    def setUp(self):
        super(InitMainFeed, self).setUp()
        self.userfeed = Stream.objects.create_main(user=self.user)

    def test_init_main_feed(self):
        user_pk = self.user.id
        res = init_main.delay(user_pk)
        self.assertTrue(res.successful())
