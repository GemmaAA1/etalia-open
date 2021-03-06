# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import datetime
from django.contrib.auth import get_user_model
from django_filters import CharFilter, MethodFilter, DateFilter, ChoiceFilter, \
    ModelMultipleChoiceFilter
from django.db.models import Q
from rest_framework import filters
from .models import Thread
from .constants import THREAD_TYPES, THREAD_PRIVACIES, THREAD_BANNED, \
    THREAD_JOINED, THREAD_LEFT, THREAD_PINNED, THREAD_INVITE_PENDING, \
    THREAD_INVITE_ACCEPTED, THREAD_PRIVATE

User = get_user_model()


class ThreadFilter(filters.FilterSet):

    doi = CharFilter(name='paper__id_doi')
    pmi = CharFilter(name='paper__id_pmi')
    arx = CharFilter(name='paper__id_arx')
    pii = CharFilter(name='paper__id_pii')
    type = ChoiceFilter(choices=THREAD_TYPES)
    privacy = ChoiceFilter(choices=THREAD_PRIVACIES)
    private = MethodFilter()
    title = CharFilter(name='title', lookup_expr='icontains')
    author = MethodFilter()
    user_id = ModelMultipleChoiceFilter(
        name='user',
        queryset=User.objects.all(),
    )
    min_date = DateFilter(name='published_at', lookup_expr='gte')
    max_date = DateFilter(name='published_at', lookup_expr='lte')
    time_span = MethodFilter(distinct=True)
    third_party = MethodFilter()

    class Meta:
        model = Thread
        fields = [
            'doi',
            'pmi',
            'arx',
            'pii',
            'type',
            'privacy',
            'title',
            'author',
            'min_date',
            'max_date',
            'time_span',
        ]
        order_by = ['published_at']

    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(ThreadFilter, self).__init__(*args, **kwargs)

    def filter_private(self, queryset, value):
        query = Q(privacy=THREAD_PRIVATE)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    # def filter_type(self, queryset, value):
    #     if value == THREAD_TYP:
    #         return queryset.filter(query)
    #     elif value in ['0',  'false']:
    #         return queryset.filter(~query)

    def filter_author(self, queryset, value):
        return queryset.filter(Q(user__first_name__icontains=value) |
                               Q(user__last_name__icontains=value))

    def filter_time_span(self, queryset, value):
        value = int(value)
        date = (datetime.datetime.now() - datetime.timedelta(days=value)).date()
        return queryset.filter(published_at__gte=date)

    def filter_third_party(self, queryset, value):
        if value in ['0',  'false']:
            third_party_models = Thread.THIRD_PARTY_TYPES
            kwargs = dict((('{0}__isnull'.format(n.lower()), True)
                           for _, n in third_party_models))
            return queryset.filter(**kwargs)
        else:
            kwargs = {'{0}__isnull'.format(value): False}
            return queryset.filter(**kwargs)


class MyThreadFilter(ThreadFilter):

    owned = MethodFilter()
    scored = MethodFilter(distinct=True)
    pinned = MethodFilter()
    banned = MethodFilter()
    joined = MethodFilter()
    left = MethodFilter()
    published = MethodFilter()
    invited = MethodFilter()
    invited_pending = MethodFilter()
    invited_accepted = MethodFilter()

    class Meta:
        model = Thread
        fields = [
            'doi',
            'pmi',
            'arx',
            'pii',
            'type',
            'privacy',
            'title',
            'author',
            'min_date',
            'max_date',
            'time_span',
            'owned',
            'scored',
            'pinned',
            'banned',
            'joined',
            'left',
            'published',
            'invited',
            'invited_pending',
            'invited_accepted',
        ]
        
    @property    
    def qs(self):
        """To deal with ordering"""
        qs = super(MyThreadFilter, self).qs
        scored = self.data.get('scored', '0')
        pinned = self.data.get('pinned', '0')
        joined = self.data.get('joined', '0')
        if scored in ['1', 'true']:
            qs = qs.order_by('-threadfeedthreads__score',
                             '-published_at')
        elif joined in ['1', 'true']:  # in my-library | papers
            qs = qs.order_by('-threaduser__modified')
        elif pinned in ['1', 'true']:  # in my-library | pinned
            qs = qs.order_by('-threaduser__modified')
        
        self._qs = qs
        return self._qs
        
    def filter_owned(self, queryset, value):
        if value in ['1',  'true']:
            return queryset.filter(user=self.request.user)
        elif value in ['0',  'false']:
            return queryset.filter(~Q(user=self.request.user))
        
    def filter_pinned(self, queryset, value):
        query = Q(threaduser__user=self.request.user) & \
                Q(threaduser__watch=THREAD_PINNED)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_banned(self, queryset, value):
        query = Q(threaduser__user=self.request.user) & \
                Q(threaduser__watch=THREAD_BANNED)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_joined(self, queryset, value):
        query = Q(threaduser__user=self.request.user) & \
                Q(threaduser__participate=THREAD_JOINED)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_left(self, queryset, value):
        query = Q(threaduser__user=self.request.user) & \
                Q(threaduser__participate=THREAD_LEFT)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_published(self, queryset, value):
        if value in ['1',  'true']:
            return queryset.filter(
                Q(user=self.request.user) &
                ~Q(published_at=None)
            )
        elif value in ['0',  'false']:
            return queryset.filter(
                Q(user=self.request.user) &
                Q(published_at=None)
            )

    def filter_invited(self, queryset, value):
        query = Q(threaduserinvite__to_user=self.request.user)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_invited_pending(self, queryset, value):
        query = Q(threaduserinvite__to_user=self.request.user) & \
                Q(threaduserinvite__status=THREAD_INVITE_PENDING)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_invited_accepted(self, queryset, value):
        query = Q(threaduserinvite__to_user=self.request.user) & \
                Q(threaduserinvite__status=THREAD_INVITE_ACCEPTED)
        if value in ['1',  'true']:
            return queryset.filter(query)
        elif value in ['0',  'false']:
            return queryset.filter(~query)

    def filter_scored(self, queryset, value):
        if value in ['1',  'true']:
            feed_name = self.data.get('feed', 'main')
            return queryset.filter(
                Q(threadfeedthreads__threadfeed__name=feed_name) &
                Q(threadfeedthreads__threadfeed__user=self.request.user)
            )
        return queryset
