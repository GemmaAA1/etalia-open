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
from users.models import UserLibPaper
from nlp.models import PaperVectors, JournalVectors, Model


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
        target_papers_pk = list(Paper.objects.filter(
            Q(date_ep__gt=from_date) |
            (Q(date_pp__gt=from_date) & Q(date_ep=None))).exclude(
            pk__in=paper_exclude_pks).values_list('pk', flat='True'))

        # computing scores
        objs_list = []
        if target_papers_pk:
            scores, pks = self.score_multi_paper(target_papers_pk)
            ind = np.argsort(scores)[::-1]
            nb_paper_max = 100
            count = 0
            while count < nb_paper_max or count == len(ind):
                objs_list.append(UserFeedPaper(
                    feed=self,
                    paper_id=pks[ind[count]],
                    score=scores[ind[count]],
                    is_score_computed=True))
                count += 1
            UserFeedPaper.objects.bulk_create(objs_list)

    def score_multi_paper(self, target_papers_pk):
        """Score a multitude of UserFeedPaper relationship (objs) at once
        """
        model_pk = self.user.settings.model.pk
        # scoring parameters
        alpha = .25
        # get vector size without padding none values as stored in db
        vector_size = Model.objects.get(id=model_pk).size

        # get seed paper pk
        seed_papers_pk = self.papers_seed.all().values('pk')

        # get related seed data
        papers_s = PaperVectors.objects\
            .filter(
                paper__pk__in=seed_papers_pk,
                model__pk=model_pk)\
            .values('vector',
                    'paper__pk',
                    'paper__date_ep',
                    'paper__date_pp',
                    'paper__journal__pk')
        jpk = list(set([sd['paper__journal__pk'] for sd in papers_s]))
        journals_s = JournalVectors.objects.filter(
            journal__pk__in=jpk,
            model__pk=model_pk).values('vector', 'journal__pk')
        journals_s_dict = dict([(jd['journal__pk'], jd['vector'][:vector_size])
                                for jd in journals_s])

        # build seed matrix
        seed_mat = np.zeros((len(papers_s), vector_size), dtype=np.float)
        for i, entry in enumerate(papers_s):
            if journals_s_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector'][:vector_size]
                journal_vec = journals_s_dict.get(entry['paper__journal__pk'])
                seed_mat[i] = (1.-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                seed_mat[i] = np.array(entry['vector'][:vector_size])

        # Build target data
        # get related target data
        papers_t = PaperVectors.objects\
            .filter(
                paper__pk__in=target_papers_pk,
                model__pk=model_pk)\
            .exclude(
            # empty absract because has unreliable signature
                paper__abstract='')\
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

        # build seed matrix
        target_mat = np.zeros((len(papers_t), vector_size),
                            dtype=np.float)
        for i, entry in enumerate(papers_t):
            if journals_t_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector'][:vector_size]
                journal_vec = journals_t_dict.get(entry['paper__journal__pk'])
                target_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                target_mat[i] = np.array(entry['vector'][:vector_size])

        pks_t = [pd['paper__pk'] for pd in papers_t]
        pks_s = [pd['paper__pk'] for pd in papers_s]

        # method 1: average of cosine similarity of journal weighted vectors)
        # scores = 1. - np.average(distance.cdist(seed_mat, target_mat, 'cosine'),
        #                          axis=0)
        # method 2: only threhold docs are counted in score
        # threshold = 0.6
        # dis = 1. - distance.cdist(seed_mat, target_mat, 'cosine')
        # dis = np.where(dis > threshold, 1, 0)
        # scores = np.sum(dis, axis=0)

        # method 3: average cosine similarity with date + journal weighted vectors
        dis = - distance.cdist(seed_mat, target_mat, 'cosine') + 1.0
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

        # Compute corresponding weights and build array
        date_wei_dic = dict([(obj['paper__pk'],
                              self.weight_date(date_days_lapse[i]))
                             for i, obj in enumerate(objs)])
        date_weights = np.zeros((len(seed_papers_pk), ))
        for i, pk in enumerate(pks_s):
            date_weights[i] = date_wei_dic[pk]

        scores = np.average(dis, weights=date_weights, axis=0)

        return scores, pks_t

    @staticmethod
    def weight_date(day_lapse):
        """Logistic function with constant term
        """
        baseline = 0.5
        k = 0.1
        delay = - 2. * 30.  # logistic function kicks off at that time from now

        return baseline + (1-baseline) / (1 + np.exp(- (day_lapse - delay) * k))

    def truc(self):
        return print('Truc\n')

    def score_paper(self, obj):
        """Score a UserFeedPaper relationship (obj)
        """
        model_pk = self.user.settings.model.pk
        # Scoring parameters
        alpha = .20
        # get vector size without padding none values as stored in db
        vector_size = Model.objects.get(id=model_pk).size

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
        journal_dict = dict([(jd['journal__pk'], jd['vector'][:vector_size])
                             for jd in journal_data])

        # build seed matrix
        seed_mat = np.zeros((len(paper_data), vector_size), dtype=np.float)
        for i, entry in enumerate(paper_data):
            if journal_dict.get(entry['paper__journal__pk'], None):  # has journal reference
                paper_vec = entry['vector'][:vector_size]
                journal_vec = journal_dict.get(entry['paper__journal__pk'])
                seed_mat[i] = (1-alpha) * np.array(paper_vec) + \
                              alpha * np.array(journal_vec)
            else:
                seed_mat[i] = np.array(entry['vector'][:vector_size])

        # Reshape target paper data
        try:
            entry = (PaperVectors.objects.get(
                         model__pk=model_pk,
                         paper__pk=obj.paper.pk).vector[:vector_size],
                     JournalVectors.objects.get(
                         model__pk=model_pk,
                         journal__pk=obj.paper.journal.pk).vector[:vector_size])
        except JournalVectors.DoesNotExist:
            entry = (PaperVectors.objects.get(
                         model__pk=model_pk,
                         paper__pk=obj.paper.pk).vector[:vector_size],
                     None)
        if entry[1]:
            target_data = np.array((1-alpha) * np.array(entry[0]) + alpha * np.array(entry[1]))
        else:
            target_data = np.array(entry[0])

        # Method #1:
        obj.score = 1. - np.average(distance.cdist(seed_mat, target_data[None, :], 'cosine'))
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


class UserFeedVector(TimeStampedModel):

    feed = models.ForeignKey(UserFeed)

    model = models.ForeignKey(Model)


