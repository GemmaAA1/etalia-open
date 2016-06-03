# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json

from django.views.generic.base import ContextMixin

from .forms import CreateUserFeedForm
from .models import StreamPapers, TrendPapers


class CreateFeedModalMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(CreateFeedModalMixin, self).get_context_data(**kwargs)
        if not self.request.user.is_anonymous():
            context['form_create_feed'] = CreateUserFeedForm()
        return context


class NavFlapMixin(object):
    """To generate context for navigation flap"""

    def get_context_data(self, **kwargs):
        context = super(NavFlapMixin, self).get_context_data(**kwargs)
        context.update(self.get_context_counters_since_last_seen())
        return context

    def get_context_counters_since_last_seen(self):

        return {
            'stream_counter': StreamPapers.objects.filter(
                stream__user=self.request.user,
                stream__name='main',
                new=True).count(),
            'trend_counter': TrendPapers.objects.filter(
                trend__user=self.request.user,
                trend__name='main',
                new=True).count(),
            'Threads_counter': 0
        }