from django.db import models
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib.auth.models import BaseUserManager
from django.db.models import Q, F
from django.utils import timezone

from core.models import TimeStampedModel
from .validators import validate_feed_name
from library.models import Paper


class UserFeedManager(BaseUserManager):
    def init_userfeed(self, name, user, papers_seed, **kwargs):
        user_feed = self.model(user=user, name=name, **kwargs)
        user_feed.save(using=self._db)
        for paper in papers_seed:
            user_feed.papers_seed.add(paper)
        user_feed.save(using=self._db)
        return user_feed

    def init_default_userfeed(self, user, **kwargs):
        """Populate a userfeed 'main' with all papers in user library
        """
        papers_seed = user.lib.papers.all()
        return self.init_userfeed('main', user, papers_seed, **kwargs)


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

    def __str__(self):
        return self.name

    def initialize(self):
        # get library matrix vector for paper seed of feed
        model_pk = self.user.settings.model.pk
        seed_papers = self.papers_seed.all()
        seed_papers_vec = [paper.papervectors_set.get(model__pk=model_pk).vector for paper in seed_papers]
        seed_journals_pk = [paper.journal.pk for paper in seed_papers if paper.journal]

        # get papers to look at
        from_date = (timezone.now() - timezone.timedelta(
            days=self.user.settings.time_lapse)).date()
        target_papers = Paper.objects.filter(Q(date_ep__gt=from_date) | (Q(date_pp__gt=from_date) & Q(date_ep=None)))

        

class UserFeedPaper(TimeStampedModel):
    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    is_in_user_lib = models.BooleanField(default=False)

    is_disliked = models.BooleanField(default=False)

    is_liked = models.BooleanField(default=False)

    class Meta:
        ordering = ['-score']
