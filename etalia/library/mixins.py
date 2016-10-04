# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db.models import Prefetch, Q

from etalia.users.constants import USER_INDIVIDUAL

from etalia.threads.models import Thread
from etalia.threads.constants import THREAD_PRIVATE
from etalia.feeds.models import StreamPapers, TrendPapers

from .models import PaperUser


class PaperEagerLoadingMixin(object):

    @staticmethod
    def setup_eager_loading(queryset, **kwargs):
        """ Perform necessary eager loading of data. """
        queryset = queryset.select_related('journal')
        queryset = queryset.prefetch_related('authors')
        queryset = queryset.prefetch_related('authorpaper_set')
        queryset = queryset.prefetch_related(
            Prefetch(
                'thread_set',
                to_attr='threads',
                queryset=Thread.objects.filter(
                    ~Q(published_at=None) & ~Q(privacy=THREAD_PRIVATE)
                )
            )
        )
        if 'user'in kwargs and \
            kwargs['user'].is_authenticated() and \
            kwargs['user'].type == USER_INDIVIDUAL:
            queryset = queryset.prefetch_related(
                Prefetch(
                    'paperuser_set',
                    to_attr='pu',
                    queryset=PaperUser.objects.filter(user=kwargs['user'])
                )
            )
            if 'request' in kwargs:
                request = kwargs.get('request')
                scored = request.query_params.get('scored')
                type = request.query_params.get('type', 'stream')
                feed_name = request.query_params.get('feed', 'main')
                if scored == '1':
                    if type == 'stream':
                        queryset = queryset.prefetch_related(
                            Prefetch(
                                'streampapers_set',
                                to_attr='sp',
                                queryset=StreamPapers.objects.filter(
                                    stream__user=kwargs['user'],
                                    stream__name=feed_name)
                            )
                        )
                    if type == 'trend':
                        queryset = queryset.prefetch_related(
                            Prefetch(
                                'trendpapers_set',
                                to_attr='tp',
                                queryset=TrendPapers.objects.filter(
                                    trend__user=kwargs['user'],
                                    trend__name=feed_name)
                            )
                        )

        if kwargs.get('altmetric'):
            queryset = queryset.filter(altmetric__isnull=False)
            queryset = queryset.select_related('altmetric__score')

        return queryset
