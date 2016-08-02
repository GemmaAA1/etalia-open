# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import django_filters
from rest_framework import filters
from ..models import Thread


class ThreadFilter(filters.FilterSet):

    doi = django_filters.CharFilter(name='paper__id_doi')
    title = django_filters.CharFilter(name='title', lookup_expr='icontains')
    min_date = django_filters.DateFilter(name='published_at', lookup_expr='gte')
    max_date = django_filters.DateFilter(name='published_at', lookup_expr='lte')

    class Meta:
        model = Thread
        fields = ['title',
                  'doi',
                  'min_date',
                  'max_date']