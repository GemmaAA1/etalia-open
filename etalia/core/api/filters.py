# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import filters


class DisabledHTMLContextualFilterBackend(filters.DjangoFilterBackend):
    """Disable HTML filter backend and pass request to Filter __init__"""

    def to_html(self, request, queryset, view):
        return ""

    def filter_queryset(self, request, queryset, view):
        filter_class = self.get_filter_class(view, queryset)

        if filter_class:
            return filter_class(request.query_params,
                                queryset=queryset,
                                request=view.request).qs

        return queryset


class DisabledHTMLSearchFilterBackend(filters.SearchFilter):

    def to_html(self, request, queryset, view):
        return ""