# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json

from django.views.generic.base import ContextMixin

from .forms import CreateUserFeedForm


class CreateFeedModalMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(CreateFeedModalMixin, self).get_context_data(**kwargs)
        if not self.request.user.is_anonymous():
            context['form_create_feed'] = CreateUserFeedForm()
        return context


class XMLFeedMixin(ContextMixin):
    """Mixin to add to context json data corresponding to controls.
    Requires BasePaperListView.
    """

    def get_context_data(self, **kwargs):
        context = super(XMLFeedMixin, self).get_context_data(**kwargs)
        context.update(self.get_context_filter_json(context))
        return context

    def get_context_filter_json(self, context):
        # build JSON object
        filter_ = {
            'filters': [],
            'pin': self.like_flag,
            'time_span': self.time_span,
            'cluster': None,
            'search_query': self.search_query,
            'number_of_papers': context['number_of_papers'],
        }

        # add to journal filter context if journal is checked
        entries = []
        for j in context['journals']:
            is_checked = j[0] in self.journals_filter
            entries.append({
                'pk': j[0],
                'name': j[1],
                'count': j[2],
                'is_checked': is_checked
            })
        filter_['filters'].append({'id': 'journal', 'entries': entries})

        # add to author filter context if author is checked
        entries = []
        for a in context['authors']:
            is_checked = a[0] in self.authors_filter
            entries.append({
                'pk': a[0],
                'name': a[1],
                'count': None,
                'is_checked': is_checked
            })
        filter_['filters'].append({'id': 'author', 'entries': entries})

        return {'filter': json.dumps(filter_)}