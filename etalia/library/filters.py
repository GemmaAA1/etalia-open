# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import datetime
from dateutil.parser import parse

from django_filters import CharFilter, MethodFilter, ModelMultipleChoiceFilter
from django.db.models import Q
from django.db.models.expressions import RawSQL
from django.conf import settings
from rest_framework import filters

from .constants import PAPER_ADDED, PAPER_TRASHED, PAPER_PINNED, PAPER_BANNED
from .models import Paper, Journal, Author


class PaperFilter(filters.FilterSet):

    doi = CharFilter(name='id_doi')
    pmi = CharFilter(name='id_pmi')
    arxiv = CharFilter(name='id_arx')
    pii = CharFilter(name='id_pii')
    title = CharFilter(name='title', lookup_expr='icontains')
    journal = MethodFilter()
    issn = MethodFilter()
    min_date = MethodFilter()
    max_date = MethodFilter()
    time_span = MethodFilter()

    class Meta:
        model = Paper
        fields = [
            'doi',
            'pmi',
            'arxiv',
            'pii',
            'title',
            'journal',
            'issn',
            'min_date',
            'max_date',
            'time_span',
        ]
        order_by = (
            ('date_fs', 'Date first seen'),
            ('altmetric__score', 'Altmetric Score')
        )

    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(PaperFilter, self).__init__(*args, **kwargs)

    def filter_journal(self, queryset, value):
        return queryset.filter(
            Q(journal__title__icontaines=value) |
            Q(journal__short_title__icontaines=value)
        )

    def filter_min_date(self, queryset, value):
        date = parse(value).date()
        return queryset.annotate(date_fp=RawSQL(
            "LEAST(date_ep, date_pp, date_fs)",
            ())
        ).filter(date_fp__gte=date)

    def filter_max_date(self, queryset, value):
        date = parse(value).date()
        return queryset.annotate(date_fp=RawSQL(
            "LEAST(date_ep, date_pp, date_fs)",
            ())
        ).filter(date_fp__lte=date)

    def filter_issn(self, queryset, value):
        return queryset.filter(
            Q(journal__id_issn=value) |
            Q(journal__id_eissn=value)
        )

    def filter_time_span(self, queryset, value):
        value = int(value)
        date = (datetime.datetime.now() - datetime.timedelta(days=value)).date()
        return queryset.annotate(date_fp=RawSQL(
            "LEAST(date_ep, date_pp, date_fs)",
            ())
        ).filter(date_fp__gte=date)


class MyPaperFilter(PaperFilter):

    added = MethodFilter()
    trashed = MethodFilter()
    pinned = MethodFilter()
    banned = MethodFilter()
    journal_id = ModelMultipleChoiceFilter(
        name='journal',
        queryset=Journal.objects.all(),
    )
    author_id = ModelMultipleChoiceFilter(
        name='authors',
        queryset=Author.objects.all(),
    )
    scored = MethodFilter()
    feed = MethodFilter()

    class Meta:
        model = Paper
        fields = [
            'doi',
            'pmi',
            'arxiv',
            'pii',
            'title',
            'journal',
            'issn',
            'min_date',
            'max_date',
            'time_span',
            'added',
            'trashed',
            'pinned',
            'banned',
            'journal_id',
            'author_id',
            'scored',
            'feed',
        ]

    def __init__(self, *args, **kwargs):
        super(MyPaperFilter, self).__init__(*args, **kwargs)
        self.form.fields['journal_id'].always_filter = False
        self.form.fields['author_id'].always_filter = False

    @property
    def qs(self):
        """To deal with ordering"""
        qs = super(MyPaperFilter, self).qs

        # Manage ordering
        type = self.data.get('type', 'stream')
        scored = self.data.get('scored', '0')
        pinned = self.data.get('pinned', '0')
        added = self.data.get('added', '0')
        if scored in ['1', 'true']:
            if type == 'stream':
                qs = qs.order_by('-streampapers__score')
            elif type == 'trend':
                qs = qs.order_by('-trendpapers__score')
        elif added in ['1', 'true']:  # in my-library | papers
            qs = qs.order_by('-userlib_paper__date_created')
        elif pinned in ['1', 'true']:  # in my-library | pinned
            qs = qs.order_by('-paperuser__modified')

        self._qs = qs
        return self._qs

    def filter_added(self, queryset, value):
        query = Q(paperuser__user=self.request.user) & \
                Q(paperuser__store=PAPER_ADDED)
        if value in ['1', 'true']:
            return queryset.filter(query)
        elif value in ['0', 'false']:
            return queryset.filter(~query)

    def filter_trashed(self, queryset, value):
        query = Q(paperuser__user=self.request.user) & \
                Q(paperuser__store=PAPER_TRASHED)
        if value in ['1', 'true']:
            return queryset.filter(query)
        elif value in ['0', 'false']:
            return queryset.filter(~query)

    def filter_pinned(self, queryset, value):
        query = Q(paperuser__user=self.request.user) & \
                Q(paperuser__watch=PAPER_PINNED)
        if value in ['1', 'true']:
            return queryset.filter(query)
        elif value in ['0', 'false']:
            return queryset.filter(~query)

    def filter_banned(self, queryset, value):
        query = Q(paperuser__user=self.request.user) & \
                Q(paperuser__watch=PAPER_BANNED)
        if value in ['1', 'true']:
            return queryset.filter(query)
        elif value in ['0', 'false']:
            return queryset.filter(~query)

    def filter_scored(self, queryset, value):
        if value in ['1', 'true']:
            type = self.data.get('type', 'stream')
            feed_name = self.data.get('feed', 'main')

            if type == 'stream':
                return queryset.filter(
                    Q(streampapers__stream__name=feed_name) &
                    Q(streampapers__stream__user=self.request.user)
                )
            elif type == 'trend':
                return queryset.filter(
                    Q(trendpapers__trend__name=feed_name) &
                    Q(trendpapers__trend__user=self.request.user)
                )

        return queryset

    def filter_time_span(self, queryset, value):
        time_span = int(value) or settings.FEED_TIME_SPAN_DEFAULT
        cutoff_datetime = (datetime.datetime.now() -
                           datetime.timedelta(days=int(time_span))).date()
        scored = bool(int(self.data.get('scored', '0')))
        type = self.data.get('type', 'stream')

        if scored:
            if type == 'stream':
                return queryset.filter(
                    Q(streampapers__date__gt=cutoff_datetime))
            elif type == 'trend':
                return queryset.filter(
                    Q(trendpapers__date__gt=cutoff_datetime))
        else:
            return queryset.filter(
                Q(userlib_paper__date_created__gt=cutoff_datetime))
