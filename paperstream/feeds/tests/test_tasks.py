from ..tasks import update_feed
from ..models import UserFeed

from .base import UserFeedTestCase

class UpdateTest(UserFeedTestCase):

    def setUp(self):
        super(UpdateTest, self).setUp()
        self.userfeed = UserFeed.objects.create_default(user=self.user)

    def test_update_async(self):
        res = update_feed.delay(self.userfeed.pk)
        self.assertTrue(res.successful())
