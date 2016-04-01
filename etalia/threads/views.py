# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic import TemplateView, DetailView, ListView
from braces.views import LoginRequiredMixin

from etalia.core.mixins import AjaxableResponseMixin
from .models import Thread
from .api.serializers import ThreadNestedSerializer


class ThreadView(LoginRequiredMixin, AjaxableResponseMixin, DetailView):

    model = Thread
    template_name = 'threads/thread.html'

    def get_context_data(self, **kwargs):
        context = super(ThreadView, self).get_context_data(**kwargs)
        context['thread'] = ThreadNestedSerializer(
            instance=self.get_object(),
            context={'request': self.request}).data
        return context

thread = ThreadView.as_view()


class MyThreadsView(LoginRequiredMixin, AjaxableResponseMixin, ListView):

    template_name = 'threads/my_threads.html'

    def get_queryset(self):
        return Thread.objects.all()

my_threads = MyThreadsView.as_view()