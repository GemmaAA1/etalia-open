# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re

from django.db import models

from altmetric import Altmetric
from altmetric import AltmetricException, AltmetricHTTPException

from django.conf import settings

from etalia.library.models import Paper
from etalia.core.models import TimeStampedModel, NullableCharField


class AltmetricModel(TimeStampedModel):

    paper = models.OneToOneField(Paper, related_name='altmetric')

    score = models.FloatField(default=0.0, db_index=True)

    # altmetric ids
    altmetric_id = models.IntegerField(default=0)
    altmetric_jid = models.CharField(max_length=32, blank=True)

    # history
    score_1d = models.FloatField(default=0.0)
    score_2d = models.FloatField(default=0.0)
    score_3d = models.FloatField(default=0.0)
    score_4d = models.FloatField(default=0.0)
    score_5d = models.FloatField(default=0.0)
    score_6d = models.FloatField(default=0.0)
    score_1w = models.FloatField(default=0.0)
    score_1m = models.FloatField(default=0.0)
    score_3m = models.FloatField(default=0.0)
    score_6m = models.FloatField(default=0.0)
    score_1y = models.FloatField(default=0.0)

    # metric details
    cited_by_posts_count = models.IntegerField(default=0)
    cited_by_delicious_count = models.IntegerField(default=0)
    cited_by_fbwalls_count = models.IntegerField(default=0)
    cited_by_feeds_count = models.IntegerField(default=0)
    cited_by_forum_count = models.IntegerField(default=0)
    cited_by_gplus_count = models.IntegerField(default=0)
    cited_by_linkedin_count = models.IntegerField(default=0)
    cited_by_msm_count = models.IntegerField(default=0)
    cited_by_peer_review_sites_count = models.IntegerField(default=0)
    cited_by_pinners_count = models.IntegerField(default=0)
    cited_by_policies_count = models.IntegerField(default=0)
    cited_by_qs_count = models.IntegerField(default=0)
    cited_by_rdts_count = models.IntegerField(default=0)
    cited_by_rh_count = models.IntegerField(default=0)
    cited_by_tweeters_count = models.IntegerField(default=0)
    cited_by_videos_count = models.IntegerField(default=0)
    cited_by_weibo_count = models.IntegerField(default=0)
    cited_by_wikipedia_count = models.IntegerField(default=0)

    readers_citeulike = models.IntegerField(default=0)
    readers_mendeley = models.IntegerField(default=0)

    # images
    image = models.URLField(default='')
    type = models.CharField(default='zzzzzzzz', max_length=8)

    class Meta:
        ordering = ['-score']

    def __str__(self):
        return '{short_title}: {score}'.format(
            short_title=self.paper.short_title,
            score=self.score)

    @property
    def print_history(self):
        s_s = []
        for field in self._meta.fields:
            if field.name.startswith('score_'):
                s_s.append((field.name, getattr(self, field.name)))
        return s_s

    @property
    def print_cited_by(self):
        cbs = []
        for field in self._meta.fields:
            if field.name.startswith('cited_by'):
                cbs.append((field.name, getattr(self, field.name)))
        return cbs

    @property
    def print_badge(self):
        url_base = 'https://altmetric-badges.a.ssl.fastly.net/?size=128&score={score:.0f}&types={type}'
        return url_base.format(score=self.score, type=self.type)

    def update(self):
        rsp = None
        a = Altmetric(settings.ALTMETRIC_API_KEY)
        ids = self.paper.get_ids()
        # Fetch altmetric data based on paper supported id
        try:
            if 'doi' in ids:
                rsp = a.doi(ids.get('doi', ''))
            if not rsp and 'arx' in ids:
                rsp = a.arxiv(ids.get('arx'))
            if not rsp and 'pmi' in ids:
                rsp = a.pmid(ids.get('pmi'))
        except AltmetricHTTPException:
            raise

        # Parse json
        if rsp:
            self.score = rsp.get('score')
            self.altmetric_id = rsp.get('altmetric_id')
            self.altmetric_jid = rsp.get('altmetric_jid', '')
            self.score_1d = rsp.get('history').get('1d')
            self.score_2d = rsp.get('history').get('2d')
            self.score_3d = rsp.get('history').get('3d')
            self.score_4d = rsp.get('history').get('4d')
            self.score_5d = rsp.get('history').get('5d')
            self.score_6d = rsp.get('history').get('6d')
            self.score_1w = rsp.get('history').get('1w')
            self.score_1m = rsp.get('history').get('1m')
            self.score_3m = rsp.get('history').get('3m')
            self.score_6m = rsp.get('history').get('6m')
            self.score_1y = rsp.get('history').get('1y')
            self.cited_by_posts_count = rsp.get('cited_by_posts_count', 0)
            self.cited_by_delicious_count = rsp.get('cited_by_delicious_count', 0)
            self.cited_by_fbwalls_count = rsp.get('cited_by_fbwalss_count', 0)
            self.cited_by_feeds_count = rsp.get('cited_by_feeds_count', 0)
            self.cited_by_forum_count = rsp.get('cited_by_forum_count', 0)
            self.cited_by_gplus_count = rsp.get('cited_by_gplus_count', 0)
            self.cited_by_linkedin_count = rsp.get('cited_by_linkedin_count', 0)
            self.cited_by_msm_count = rsp.get('cited_by_msm_count', 0)
            self.cited_by_peer_review_sites_count = rsp.get('cited_by_peer_review_sites_count', 0)
            self.cited_by_pinners_count = rsp.get('cited_by_pinners_count', 0)
            self.cited_by_policies_count = rsp.get('cited_by_policies_count', 0)
            self.cited_by_qs_count = rsp.get('cited_by_qs_count', 0)
            self.cited_by_rdts_count = rsp.get('cited_by_rdts_count', 0)
            self.cited_by_rh_count = rsp.get('cited_by_rh_count', 0)
            self.cited_by_tweeters_count = rsp.get('cited_by_tweeters_count', 0)
            self.cited_by_videos_count = rsp.get('cited_by_videos_count', 0)
            self.cited_by_weibo_count = rsp.get('cited_by_weibo_count', 0)
            self.cited_by_wikipedia_count = rsp.get('cited_by_wikipedia_count', 0)
            if isinstance(rsp.get('readers'), dict):
                self.readers_citeulike = rsp.get('readers').get('citeulike')
                self.readers_mendeley = rsp.get('readers').get('mendeley')
            self.image = rsp.get('images').get('small')
            match = re.findall(r'types=(P?[\w]+)', self.image)
            if match:
                self.type = match[0]

        # save
        self.save()






