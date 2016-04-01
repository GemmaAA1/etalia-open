# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic import UpdateView, FormView, TemplateView, CreateView, \
    DeleteView, ListView
from braces.views import LoginRequiredMixin, UserPassesTestMixin
from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404

from etalia.core.mixins import AjaxableResponseMixin
from .forms import ThreadCreateForm, ThreadUpdateForm, ThreadPostForm, \
    ThreadPostCommentForm, ThreadUserForm
from .models import Thread, ThreadPost, ThreadComment, ThreadUser
from .api.serializers import ThreadSerializer, ThreadPostSerializer, \
    ThreadCommentSerializer, ThreadUserSerializer


class ThreadView(LoginRequiredMixin, AjaxableResponseMixin, TemplateView):

    template_name = 'threads/thread.html'

thread = ThreadView.as_view()


class MyThreadsView(LoginRequiredMixin, AjaxableResponseMixin, TemplateView):

    template_name = 'threads/my_threads.html'

    def get_queryset(self):
        return Thread.objects.all()

my_threads = MyThreadsView.as_view()