# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db.models import Prefetch, Q

from etalia.users.constants import USER_INDIVIDUAL

from .models import PaperUser
from etalia.threads.models import Thread
from etalia.threads.constant import THREAD_PRIVATE


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
        return queryset
