from operator import add

import numpy as np
from scipy.spatial import distance

from django.db import models
from django.contrib.auth import get_user_model
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

from config.celery import celery_app as app


class UserFeedManager(BaseUserManager):

    def create(self, **kwargs):
        """Create new UserFeed / UserFeedVector
        """
        if 'name' not in kwargs:
            raise AssertionError('<name> is not defined')
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')

        if 'papers_seed' in kwargs:
            papers_seed = kwargs['papers_seed']
            kwargs.pop('papers_seed', None)
        else:
            papers_seed = None
        user_feed = super(UserFeedManager, self).create(**kwargs)

        # Init papers_seed
        user_feed.papers_seed.add(*papers_seed)

        # Init feed vector
        model_pks = Model.objects.all().values_list('pk', flat='True')
        for model_pk in model_pks:
            ufv = UserFeedVector.objects.create(model_id=model_pk,
                                                feed=user_feed)
            ufv.update_vector()

        # update feed
        user_feed.update()

        return user_feed

    def create_default(self, **kwargs):
        """Create a new default UserFeed / UserFeedVector, named 'main' with
        all papers in user library
        """
        if 'user' not in kwargs:
            raise AssertionError('<user> is not defined')
        papers = kwargs.get('user').lib.papers.all()
        return self.create(name='main', papers_seed=papers, **kwargs)


class UserFeed(TimeStampedModel):
    """User Feed model

    Feed are based on what user has defined as papers_seed
    """

    name = models.CharField(max_length=100, default='main',
                            validators=[validate_feed_name])

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='feed')

    # Papers that are used in the similarity matching
    papers_seed = models.ManyToManyField(Paper,
                                         related_name='feed_seed')

    # relevant papers matched
    papers_match = models.ManyToManyField(Paper, through='UserFeedPaper',
                                          related_name='feed_match')

    status = models.CharField(max_length=3, blank=True, default='',
                              choices=(('', 'Uninitialized'),
                                       ('IDL', 'Idle'),
                                       ('ING', 'Syncing')))

    objects = UserFeedManager()

    @property
    def count_papers_seed(self):
        return self.papers_seed.all().count()

    @property
    def count_papers_match(self):
        return self.papers_match.all().count()

    def set_feed_syncing(self):
        self.feed_status = 'ING'
        self.save()

    def set_feed_idle(self):
        self.status = 'IDL'
        self.save()

    def __str__(self):
        return '{feed}@{username}'.format(feed=self.name, username=self.user.email)

    def update(self):
        """
        - Get all papers in time range excluding untrusted and already in user lib
        #TODO: Restricted this query to the k-NN to the feed vector
        - Score papers
        - Get only the N top scored papers
        - Create/Update UserFeedPaper
        """
        # time range
        from_date = (timezone.now() -
                     timezone.timedelta(
                         days=self.user.settings.time_lapse)).date()

        # Clean old papers
        UserFeedPaper.objects\
            .filter(Q(feed=self) &
                    (Q(paper__date_ep__lt=from_date) |
                    (Q(paper__date_pp__lt=from_date) & Q(paper__date_ep=None))))\
            .delete()
        ufp_pks_to_update = UserFeedPaper.objects\
            .filter(feed=self)\
            .values_list('pk', flat='True')

        # get target papers excluding user lib + not trusted + empty abstract
        paper_exclude_pks = list(self.user.lib.papers.values_list('pk',
                                                                  flat='True'))

        # Get nearest neighbors of seed papers and retrieve target_paper
        lsh_pk = LSH.objects.get(model_id=self.user.settings.model.pk,
                                 time_lapse=self.user.settings.time_lapse)
        target_papers_pk = PaperNeighbors.objects\
            .filter(paper__pk__in=self.papers_seed.values('pk'),
                    lsh_id=lsh_pk)\
            .exclude(Q(pk__in=paper_exclude_pks) |
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

    def score_multi_papers(self, target_papers_pk, scoring_method=1):
        """Score papers
        """

        model_pk = self.user.settings.model.pk
        # scoring parameters
        alpha = settings.FEED_JOURNAL_VECTOR_RATIO
        # get true vector size
        vector_size = Model.objects.get(id=model_pk).size

        # get seed paper pk
        seed_papers_pk = self.papers_seed.values('pk')

        # get seed data
        papers_s = PaperVectors.objects\
            .filter(
                paper__pk__in=seed_papers_pk,
                model__pk=model_pk)\
            .values('vector',
                    'paper__pk',
                    'paper__date_ep',
                    'paper__date_pp',
                    'paper__journal__pk')
        # get journal of seed papers
        jpk = list(set([sd['paper__journal__pk'] for sd in papers_s]))
        journals_s = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        # rearrange in a dict
        journals_s_dict = dict([(jd['journal__pk'], jd['vector'][:vector_size])
                                for jd in journals_s])

        # build seed matrix
        # init seed matrix
        seed_mat = np.zeros((len(papers_s), vector_size), dtype=np.float)
        # populate matrix
        for i, entry in enumerate(papers_s):
            # if paper has journal reference, weight paper with journal
            if journals_s_dict.get(entry['paper__journal__pk'], None):
                paper_vec = entry['vector'][:vector_size]
                journal_vec = journals_s_dict.get(entry['paper__journal__pk'])
                seed_mat[i] = (1.-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                seed_mat[i] = np.array(entry['vector'][:vector_size])

        # Build target data (same as previous with target data)
        # get related target data
        papers_t = PaperVectors.objects\
            .filter(
                paper__pk__in=target_papers_pk,
                model__pk=model_pk)\
            .values('paper__pk',
                    'vector',
                    'paper__date_ep',
                    'paper__date_pp',
                    'paper__journal__pk')
        jpk = list(set([sd['paper__journal__pk'] for sd in papers_t]))
        journals_t = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        journals_t_dict = dict([(jd['journal__pk'], jd['vector'][:vector_size])
                                for jd in journals_t])

        # build target matrix
        target_mat = np.zeros((len(papers_t), vector_size), dtype=np.float)
        for i, entry in enumerate(papers_t):
            if journals_t_dict.get(entry['paper__journal__pk'], None):
                paper_vec = entry['vector'][:vector_size]
                journal_vec = journals_t_dict.get(entry['paper__journal__pk'])
                target_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                target_mat[i] = np.array(entry['vector'][:vector_size])

        # These are the ordered pk
        pks_t = [pd['paper__pk'] for pd in papers_t]
        pks_s = [pd['paper__pk'] for pd in papers_s]

        # scoring_method #1: average of cosine similarity of journal weighted vectors)
        if scoring_method == 1:
            scores = 1. - np.average(distance.cdist(seed_mat, target_mat, 'cosine'),
                                     axis=0)
        # scoring_method #2: only threshold docs are counted in score
        elif scoring_method == 2:
            threshold = 0.6
            dis = 1. - distance.cdist(seed_mat, target_mat, 'cosine')
            dis = np.where(dis > threshold, 1, 0)
            scores = np.sum(dis, axis=0)

        # scoring_method 3: average cosine similarity with date + journal weighted
        # vectors
        elif scoring_method == 3:

            dis = 1.0 - distance.cdist(seed_mat, target_mat, 'cosine')
            # get date when paper was created in user lib
            objs = UserLibPaper.objects\
                .filter(
                    userlib=self.user.lib,
                    paper_id__in=seed_papers_pk)\
                .values(
                    'paper__pk',
                    'date_created',
                    'date_last_modified')

            # Convert date in lapse day array
            date_val = np.array([obj['date_created'] for obj in objs])
            date_days_lapse = \
                np.array(
                    list(
                        map(
                            lambda x: x.days,
                            date_val - timezone.datetime.today().date()
                        )
                    )
                )

            # Compute date weight
            date_wei_dic = dict(
                [(obj['paper__pk'], self.weight_date(date_days_lapse[i]))
                 for i, obj in enumerate(objs)])

            # Rearrange as np.array
            date_weights = np.zeros((len(seed_papers_pk), ))
            for i, pk in enumerate(pks_s):
                date_weights[i] = date_wei_dic[pk]

            # Compute score as a weighted average
            scores = np.average(dis, weights=date_weights, axis=0)

        else:
            raise ValueError('Scoring method unknown: {0}'.format(scoring_method))

        return scores, pks_t

    @staticmethod
    def weight_date(day_lapse):
        """Logistic function with constant baseline term
        """
        baseline = 0.5
        k = 0.1
        delay = - 2. * 30.  # logistic function kicks off at that time from now

        return baseline + (1-baseline) / (1 + np.exp(- (day_lapse - delay) * k))


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

