from django.db import models
from django.contrib.auth import get_user_model
# Create your models here.

from library.models import Paper

User = get_user_model()


class UserFeed(models.Model):
    """Feed of user"""

    user = models.ForeignKey(User, related_name='feed')

    papers = models.ManyToManyField(Paper, through='UserFeedPaper')

    feed_status = models.CharField(max_length=3, blank=True, default='',
                                   choices=(('', 'Uninitialized'),
                                            ('IDL', 'Idle'),
                                            ('ING', 'Syncing')))

    def count_papers(self):
        return self.papers.all().count()

    def set_feed_syncing(self):
        self.feed_status = 'ING'
        self.save()

    def set_feed_idle(self):
        self.feed_status = 'IDL'
        self.save()


class UserFeedPaper(models.Model):

    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)