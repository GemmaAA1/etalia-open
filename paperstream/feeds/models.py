from django.db import models
from django.contrib.auth import get_user_model
from django.conf import settings

from core.models import TimeStampedModel
from .validators import validate_feed_name
from library.models import Paper


class UserFeed(TimeStampedModel):
    """Feed of user"""

    name = models.CharField(max_length=100, default='Main',
                            validators=[validate_feed_name])

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='feed')

    # cluster of paper that is used in the similarity matching
    paper_seed = models.ManyToManyField(Paper, related_name='paper_in')

    # relevant papers matched
    paper_match = models.ManyToManyField(Paper, through='UserFeedPaper',
                                       related_name='paper_out')

    status = models.CharField(max_length=3, blank=True, default='',
                               choices=(('', 'Uninitialized'),
                                        ('IDL', 'Idle'),
                                        ('ING', 'Syncing')))

    @property
    def count_paper_in(self):
        return self.paper_in.all().count()

    @property
    def count_paper_out(self):
        return self.paper_out.all().count()

    def set_feed_syncing(self):
        self.feed_status = 'ING'
        self.save()

    def set_feed_idle(self):
        self.status = 'IDL'
        self.save()


class UserFeedPaper(TimeStampedModel):

    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    class Meta:
        ordering = ['-score']