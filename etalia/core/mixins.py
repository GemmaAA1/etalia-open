# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from progressbar import ProgressBar, Percentage, Bar, ETA

from django.http import JsonResponse
from django.db import models

from etalia.users.mixins import ProfileModalFormsMixin
from etalia.feeds.mixins import CreateFeedModalMixin
from etalia.feeds.models import StreamMatches, TrendMatches
from etalia.threads.models import Thread

from .models import ThreadVectors


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

    def get_ajax_data(self, *args, **kwargs):
        raise NotImplementedError


class ModalMixin(ProfileModalFormsMixin, CreateFeedModalMixin):
    """Pull Mixin in one"""
    pass


class NavFlapMixin(object):
    """To generate context for navigation flap"""

    def get_context_data(self, **kwargs):
        context = super(NavFlapMixin, self).get_context_data(**kwargs)
        context.update(self.get_context_counters_since_last_seen())
        return context

    def get_context_counters_since_last_seen(self):

        return {
            'stream_counter': StreamMatches.objects.filter(
                stream__user=self.request.user,
                stream__name='main',
                new=True).count(),
            'trend_counter': TrendMatches.objects.filter(
                trend__user=self.request.user,
                trend__name='main',
                new=True).count(),
            'Threads_counter': 0
        }


class XMLMixin(object):
    """Mixin to add to context json data corresponding to controls.

    Requires Class BasePaperListView.
    """

    def get_context_data(self, **kwargs):
        context = super(XMLMixin, self).get_context_data(**kwargs)
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


class ModelThreadMixin(object):
    """Mixin to Model to add Thread support"""
    THREAD_TOKENIZED_FIELDS = ['title', 'content']

    def infer_thread(self, thread_pk, **kwargs):
        """
        Infer vector for model and thread instance
        """
        vector = self.infer_object(Thread, thread_pk,
                                   self.THREAD_TOKENIZED_FIELDS, **kwargs)

        pv, _ = ThreadVectors.objects.get_or_create(model=self,
                                                    thread_id=thread_pk)

        # store
        pv.set_vector(vector.tolist())

        return thread_pk

    def infer_threads(self, thread_pks, **kwargs):
        """Infer vector for model and all thread in thread_pks list
        """
        # Check inputs
        if isinstance(thread_pks, models.QuerySet):
            thread_pks = list(thread_pks)
        if not isinstance(thread_pks, list):
            raise TypeError(
                'thread_pks must be list or QuerySet, found {0} instead'.format(
                    type(thread_pks)))

        # setup progressbar
        nb_pbar_updates = 100
        nb_threads = len(thread_pks)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=nb_threads, redirect_stderr=True).start()

        for count, thread_pk in enumerate(thread_pks):
            self.infer_thread(thread_pk, **kwargs)
            if not (nb_threads // nb_pbar_updates) == 0:
                if not count % (nb_threads // nb_pbar_updates):
                    pbar.update(count)
        # close progress bar
        pbar.finish()

    def tasks(self, task, **kwargs):
        if task == 'infer_thread':
            thread_pk = kwargs.pop('thread_pk')
            return self.infer_thread(thread_pk, **kwargs)
        if task == 'infer_threads':
            thread_pks = kwargs.pop('thread_pks')
            return self.infer_threads(thread_pks, **kwargs)

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented
