# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.template.response import TemplateResponse
from django.views.generic import DetailView, RedirectView
from django.utils.text import slugify
from .models import Thread


def my_threads(request):
    return TemplateResponse(
        request,
        'threads/my_list.html',
        {'control_states': json.dumps(request.session.get(
            'threads-control-states',
            {'time_span': None, 'search': None, 'pin': 0}))}
    )


def threads(request):
    return TemplateResponse(
        request,
        'threads/list.html',
        {}
    )


class ThreadDetail(DetailView):

    template_name = 'threads/my_list.html'
    model = Thread

    def get_template_names(self):
        if self.request.user.is_anonymous():
            self.template_name = 'trends/list.html'
        return super(ThreadDetail, self).get_template_names()

thread_slug = ThreadDetail.as_view()


class ThreadViewPk(RedirectView):
    """Redirect to slug paper url"""

    permanent = False
    query_string = True
    pattern_name = 'threads:thread-slug'

    def get_redirect_url(self, *args, **kwargs):
        thread = Thread.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(thread.title)
        return super(ThreadViewPk, self).get_redirect_url(*args, **kwargs)

thread = ThreadViewPk.as_view()