# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import django_filters
from rest_framework import filters
from ..models import Paper


class PaperFilter(filters.FilterSet):
    doi = django_filters.CharFilter(name='id_doi')
    journal = django_filters.CharFilter(name='journal__title', lookup_expr='icontains')

    class Meta:
        model = Paper
        fields = ['title', 'journal', 'doi']