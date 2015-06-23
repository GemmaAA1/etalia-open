from django.db import models
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib.auth.models import BaseUserManager

from core.models import TimeStampedModel
from .validators import validate_feed_name
from library.models import Paper


class UserFeedManager(BaseUserManager):

    def init_userfeed(self, name, user, papers_seed):
        uf = UserFeed(user=user, name=name)
        for paper in papers_seed:
            uf.papers_seed.add(paper)
        uf.save()
        return uf

    def init_default_userfeed(self, user):
        return self.init_userfeed('main', user, user.lib.papers.all())


class UserFeed(TimeStampedModel):
    """Feed of user"""

    name = models.CharField(max_length=100, default='main',
                            validators=[validate_feed_name])

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='feed')

    # cluster of paper that is used in the similarity matching
    papers_seed = models.ManyToManyField(Paper, related_name='paper_in')

    # relevant papers matched
    papers_match = models.ManyToManyField(Paper, through='UserFeedPaper',
                                       related_name='paper_out')

    status = models.CharField(max_length=3, blank=True, default='',
                               choices=(('', 'Uninitialized'),
                                        ('IDL', 'Idle'),
                                        ('ING', 'Syncing')))

    objects = UserFeedManager()

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

    is_in_user_lib = models.BooleanField(default=False)

    is_disliked = models.BooleanField(default=False)

    is_liked = models.BooleanField(default=False)

    class Meta:
        ordering = ['-score']