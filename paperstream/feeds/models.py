import numpy as np
from scipy.spatial import distance

from django.db import models
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib.auth.models import BaseUserManager
from django.db.models import Q, F
from django.utils import timezone

from core.models import TimeStampedModel
from .validators import validate_feed_name
from library.models import Paper
from nlp.models import PaperVectors, JournalVectors


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
        return self.papers_seed.all().count()

    @property
    def count_paper_out(self):
        return self.papers_match.all().count()

    def set_feed_syncing(self):
        self.feed_status = 'ING'
        self.save()

    def set_feed_idle(self):
        self.status = 'IDL'
        self.save()

    def __str__(self):
        return '{feed}@{username}'.format(feed=self.name, username=self.user.email)

    def initialize(self):
        """Initialize a news feed
        """
        # get papers to look at excluding papers already in user lib
        from_date = (timezone.now() - timezone.timedelta(days=self.user.settings.time_lapse)).date()
        paper_exclude_pks = self.user.lib.papers.values_list('pk', flat='True')
        target_papers_pk = Paper.objects.filter(
            Q(date_ep__gt=from_date) |
            (Q(date_pp__gt=from_date) & Q(date_ep=None))).exclude(
            pk__in=paper_exclude_pks).values('pk')

        # create related UserFeedPaper objects
        objs_list = [UserFeedPaper(feed=self, paper_id=pk)
                     for pk in target_papers_pk]
        objs = UserFeedPaper.objects.bulk_create(objs_list)

        # compute scores
        if objs:
            self.score_multi_paper(objs)

    def update_feed(self):
        # Init papers to look at
        from_date = (timezone.now() -
                     timezone.timedelta(
                         days=self.user.settings.time_lapse)).date()
        # delete old papers
        UserFeedPaper.objects.filter(
            Q(feed=self) &
            (Q(paper__date_ep__gt=from_date) |
             (Q(paper__date_pp__gt=from_date) & Q(paper__date_ep=None)))).delete()

        # get targeted papers excluding papers
        paper_exclude_pks = list(self.user.lib.papers.values_list('pk',
                                                                  flat='True'))
        paper_exclude_pks += list(self.papers_match.values_list('pk',
                                                                flat='True'))
        paper_exclude_pks = list(set(paper_exclude_pks))
        target_papers_pk = Paper.objects.filter(
            Q(date_ep__gt=from_date) |
            (Q(date_pp__gt=from_date) & Q(date_ep=None))).exclude(
            pk__in=paper_exclude_pks).values_list('pk', flat='True')

        # bulk creating new objects
        objs_list = [UserFeedPaper(feed=self, paper_id=pk)
                     for pk in target_papers_pk]
        UserFeedPaper.objects.bulk_create(objs_list)

        objs = UserFeedPaper.objects.filter(feed=self)

        # computing scores
        if objs:
            self.score_multi_paper(objs)

    def score_multi_paper(self, objs):
        """Score a multitude of UserFeedPaper relationship (objs) at once
        """
        model_pk = self.user.settings.model.pk
        # scoring parameters
        alpha = .20

        # get seed paper pk
        seed_papers_pk = self.papers_seed.all().values('pk')

        # get related seed data
        paper_data = PaperVectors.objects.filter(
            paper__pk__in=seed_papers_pk,
            model__pk=model_pk).values('vector',
                                       'paper__date_ep',
                                       'paper__date_pp',
                                       'paper__journal__pk')
        jpk = list(set([sd['paper__journal__pk'] for sd in paper_data]))
        journal_data = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        journal_dict = dict([(jd['journal__pk'], jd['vector']) for jd in journal_data])

        # build seed matrix
        seed_mat = np.zeros((len(paper_data), len(paper_data[0]['vector'])),
                            dtype=np.float)
        for i, entry in enumerate(paper_data):
            if journal_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector']
                journal_vec = journal_dict.get(entry['paper__journal__pk'])
                seed_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                seed_mat[i] = np.array(entry['vector'])

        # build target data
        # get target data pk
        target_papers_pk = objs.values('paper__pk')

        # get related target data
        paper_data = PaperVectors.objects.filter(
            paper__pk__in=target_papers_pk,
            model__pk=model_pk).values('vector',
                                       'paper__date_ep',
                                       'paper__date_pp',
                                       'paper__journal__pk')
        jpk = list(set([sd['paper__journal__pk'] for sd in paper_data]))
        journal_data = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        journal_dict = dict([(jd['journal__pk'], jd['vector']) for jd in journal_data])

        # build seed matrix
        target_mat = np.zeros((len(paper_data), len(paper_data[0]['vector'])),
                            dtype=np.float)
        for i, entry in enumerate(paper_data):
            if journal_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector']
                journal_vec = journal_dict.get(entry['paper__journal__pk'])
                target_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                target_mat[i] = np.array(entry['vector'])

        # compute score
        scores = np.sum(distance.cdist(seed_mat, target_mat, 'cosine'), axis=0)
        for i, obj in enumerate(objs):
            obj.score = scores[i]
            obj.is_score_computed = True
            obj.save()

    def score_paper(self, obj):
        """Score a UserFeedPaper relationship (obj)
        """
        model_pk = self.user.settings.model.pk
        # Scoring parameters
        alpha = .20

        # Get seed papers pk
        seed_papers_pk = self.papers_seed.all().values('pk')

        # Get related seed data
        paper_data = PaperVectors.objects.filter(
            paper__pk__in=seed_papers_pk,
            model__pk=model_pk).values('vector',
                                       'paper__date_ep',
                                       'paper__date_pp',
                                       'paper__journal__pk')
        jpk = list(set([sd['paper__journal__pk'] for sd in paper_data]))
        journal_data = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        journal_dict = dict([(jd['journal__pk'], jd['vector']) for jd in journal_data])

        # build seed matrix
        seed_mat = np.zeros((len(paper_data), len(paper_data[0]['vector'])),
                            dtype=np.float)
        for i, entry in enumerate(paper_data):
            if journal_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector']
                journal_vec = journal_dict.get(entry['paper__journal__pk'])
                seed_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                seed_mat[i] = np.array(entry['vector'])

        # Reshape target paper data
        try:
            entry = (PaperVectors.objects.get(
                         model__pk=model_pk,
                         paper__pk=obj.paper.pk).vector,
                     JournalVectors.objects.get(
                         model__pk=model_pk,
                         journal__pk=obj.paper.journal.pk).vector)
        except JournalVectors.DoesNotExist:
            entry = (PaperVectors.objects.get(
                         model__pk=model_pk,
                         paper__pk=obj.paper.pk).vector,
                     None)
        if entry[1]:
            target_data = np.array((1-alpha) * np.array(entry[0]) + alpha * np.array(entry[1]))
        else:
            target_data = np.array(entry[0])

        # Method #1:
        obj.score = np.sum(distance.cdist(seed_mat, target_data[None, :], 'cosine'))
        obj.is_score_computed = True
        obj.save()


class UserFeedPaper(TimeStampedModel):
    feed = models.ForeignKey(UserFeed)

    paper = models.ForeignKey(Paper)

    score = models.FloatField(default=0.)

    is_score_computed = models.BooleanField(default=False)

    is_disliked = models.BooleanField(default=False)

    is_liked = models.BooleanField(default=False)

    class Meta:
        ordering = ['-score']
