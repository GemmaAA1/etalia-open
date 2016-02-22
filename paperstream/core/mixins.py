# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.http import JsonResponse
from django.views.generic.base import ContextMixin
from paperstream.users.mixins import ProfileModalFormsMixin
from paperstream.feeds.mixins import CreateFeedModalMixin
from paperstream.feeds.models import StreamMatches, TrendMatches

class AjaxableResponseMixin(object):
    """
    Mixin to add AJAX support to a form.
    Must be used with an object-based FormView (e.g. UpdateView)
    """

    def form_invalid(self, form):
        response = super(AjaxableResponseMixin, self).form_invalid(form)
        if self.request.is_ajax():
            return JsonResponse(form.errors, status=400)
        else:
            return response

    def form_valid(self, form):
        # We make sure to call the parent's form_valid() method because
        # it might to do some processing (in the case of CreateView, it will
        # call form.save() for example)
        response = super(AjaxableResponseMixin, self).form_valid(form)
        if self.request.is_ajax():
            data = self.get_ajax_data(form=form)
            return JsonResponse(data)
        else:
            return response

    def get_ajax_data(self, *kwargs):
        raise NotImplementedError


class ModalMixin(ProfileModalFormsMixin, CreateFeedModalMixin):
    """Pull Mixin in one"""
    pass


class NavFlapMixin(ContextMixin):
    """To generate context for navigation flap"""

    def get_context_data(self, **kwargs):
        context = super(NavFlapMixin, self).get_context_data(**kwargs)
        context.update(self.get_context_counters_since_last_login())
        return context

    def get_context_counters_since_last_login(self):

        return {
            'stream_counter': StreamMatches.objects.filter(
                stream__user=self.request.user,
                stream__name='main',
                created__gt=self.request.user.last_login).count(),
            'trend_counter': TrendMatches.objects.filter(
                trend__user=self.request.user,
                trend__name='main',
                created__gt=self.request.user.last_login).count(),
            'tocles_counter': 0}


