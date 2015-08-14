import numpy as np

from django.db import models

from django.conf import settings
from django.contrib.auth.models import BaseUserManager
from django.db.models import Q, F
from django.utils import timezone
from django.contrib.postgres.fields import ArrayField

from core.models import TimeStampedModel
from core.utils import pad_vector
from library.models import Paper
from users.models import UserLibPaper
from nlp.models import PaperVectors, JournalVectors, Model, PaperNeighbors, LSH

from .validators import validate_feed_name
from .constants import FEED_STATUS_CHOICES
from .scoring import SimpleAverage, ThresholdAverage, WeightedJournalAverage, \
    WeightedJournalCreatedDateAverage


class UserFeedManager(BaseUserManager):

    def create(self, **kwargs):
        """Create new UserFeed / UserFeedVector
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')

        if 'papers_seed' in kwargs:
            papers_seed = kwargs['papers_seed']
            kwargs.pop('papers_seed', None)
        else:
            papers_seed = None
        obj = super(UserFeedManager, self).create(**kwargs)

        # Init papers_seed
        obj.add_seed_papers(papers_seed)

        # Init feed vector
        obj.update_userfeed_vector()

        # update feed
        obj.update()

        return obj

    def create_default(self, **kwargs):
        """Create a new default UserFeed / UserFeedVector, named 'main' with
        all papers in user library
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')
        papers = kwargs.get('user').lib.papers.all()
        return self.create(name='main', papers_seed=papers, **kwargs)


class UserFeed(TimeStampedModel):
    """User Feed model"""

    name = models.CharField(max_length=100, default='main',
                            validators=[validate_feed_name])

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='feed')

    # Papers that are used in the similarity matching
    papers_seed = models.ManyToManyField(Paper,
                                         related_name='feed_seed')

    # relevant papers matched
    papers_match = models.ManyToManyField(Paper, through='UserFeedPaper',
                                          related_name='feed_match')

    state = models.CharField(max_length=3, blank=True, default='NON',
                             choices=FEED_STATUS_CHOICES)

    objects = UserFeedManager()

    class Meta:
        unique_together = (('name', 'user'),)

    def set_state(self, state):
        self.state = state
        self.save()

    def __str__(self):
        return '{feed}@{username}'.format(feed=self.name, username=self.user.email)

    def add_seed_papers(self, papers):
        if papers:
            self.papers_seed.add(*papers)

    def update_userfeed_vector(self):
        model_pks = Model.objects.all().values_list('pk', flat='True')
        for model_pk in model_pks:
            ufv = UserFeedVector.objects.create(model_id=model_pk,
                                                feed=self)
            ufv.update_vector()

    def clean_old_papers(self):
        # get time lapse corresponding to user settings
        from_date = (timezone.now() - timezone.timedelta(
            days=self.user.settings.time_lapse)).date()

        # clean old papers
        UserFeedPaper.objects\
            .filter(Q(feed=self) &
                    (Q(paper__date_ep__lt=from_date) |
                    (Q(paper__date_pp__lt=from_date) & Q(paper__date_ep=None))))\
            .delete()

    def update(self):
        """Update UserFeed

        - Get all papers in time range excluding untrusted and already in user lib
        #TODO: Restricted this query to the k-NN to the feed vector
        - Score papers
        - Get only the N top scored papers
        - Create/Update UserFeedPaper
        """

        self.clean_old_papers()

        # get paper that remain and need to be updated
        ufp_pks_to_update = self.papers_match.all().values_list('pk',
                                                                flat='True')

        # get target papers (excluding userlib + not trusted + empty abstract)
        # from LSH

        #
        lsh_pk = LSH.objects.get(model_id=self.user.settings.model.pk,
                                 time_lapse=self.user.settings.time_lapse)
        target_papers_pk = PaperNeighbors.objects\
            .filter(paper__pk__in=self.papers_seed.values('pk'),
                    lsh_id=lsh_pk)\
            .exclude(Q(pk__in=self.user.lib.papers.values('pk')) |
                     Q(is_trusted=False) |
                     Q(abstract=''))\
            .values_list('neighbors', flat=True)

        # de-nest list and make unique
        target_papers_pk = [a for b in target_papers_pk for a in b]
        target_papers_pk = list(set(target_papers_pk))

        # Match target paper
        objs_list = []
        if target_papers_pk:
            # compute scores
            scores, pks = self.score_multi_papers(
                target_papers_pk,
                scoring_method=self.user.settings.scoring_method)
            # sort scores
            ind = np.argsort(scores)[::-1]
            # create/update UserFeedPaper
            for i in range(settings.FEED_SCORE_KEEP_N_PAPERS):
                # update
                if pks[ind[i]] in ufp_pks_to_update:
                    ufp = UserFeedPaper(feed=self, paper_id=pks[ind[i]])
                    ufp.score = scores[ind[i]]
                    ufp.is_score_computed = True
                    ufp.save()
                else:  # create in bulk
                    objs_list.append(UserFeedPaper(
                        feed=self,
                        paper_id=pks[ind[i]],
                        score=scores[ind[i]],
                        is_score_computed=True))
            UserFeedPaper.objects.bulk_create(objs_list)



class UserFeedPaper(TimeStampedModel):
    """Relationship table between model and paper with scoring
    """
    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    is_score_computed = models.BooleanField(default=False)

    is_disliked = models.BooleanField(default=False)

    is_liked = models.BooleanField(default=False)

    class Meta:
        ordering = ['-score']
        unique_together = [('feed', 'paper')]

class UserFeedVector(TimeStampedModel):
    """Feature vector for feed is defined as the averaged of paper vectors in
    feed
    """
    feed = models.ForeignKey(UserFeed, related_name='vectors')

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]

    def update_vector(self):
        pv_pks = self.feed.papers_seed.all().values('vectors__pk')
        vectors = list(PaperVectors.objects.filter(id__in=pv_pks)
                       .values_list('vector', flat='true'))

        # sum of vectors
        vector = [sum(x) for x in zip(*vectors)]
        self.set_vector(vector)

    class Meta:
        unique_together = [('feed', 'model'), ]

    def __str__(self):
        return '{feed_name}/{model_name}'.format(feed_name=self.feed.name,
                                                 model_name=self.model.name)

