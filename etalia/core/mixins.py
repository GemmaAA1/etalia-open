# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json

from django.http import JsonResponse
from django.forms.models import model_to_dict

# from etalia.feeds.models import StreamPapers, TrendPapers


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


class NavFlapMixin(object):
    """To generate context for navigation flap"""

    def get_context_data(self, **kwargs):
        context = super(NavFlapMixin, self).get_context_data(**kwargs)
        context.update(self.get_context_counters_since_last_seen())
        return context

    def get_context_counters_since_last_seen(self):

        return {
            'stream_counter': 0,
            'trend_counter': 0,
            'Threads_counter': 0
        }


class ModelDiffMixin(object):
    """
    A model mixin that tracks model fields' values and provide some useful api
    to know what fields have been changed.
    """

    def __init__(self, *args, **kwargs):
        super(ModelDiffMixin, self).__init__(*args, **kwargs)
        self.__initial = self._dict

    @property
    def diff(self):
        d1 = self.__initial
        d2 = self._dict
        diffs = [(k, (v, d2[k])) for k, v in d1.items() if v != d2[k]]
        return dict(diffs)

    @property
    def has_changed(self):
        return bool(self.diff)

    @property
    def changed_fields(self):
        return self.diff.keys()

    def get_field_diff(self, field_name):
        """
        Returns a diff for field if it's changed and None otherwise.
        """
        return self.diff.get(field_name, None)

    def save(self, *args, **kwargs):
        """
        Saves model and set initial state.
        """
        super(ModelDiffMixin, self).save(*args, **kwargs)
        self.__initial = self._dict

    @property
    def _dict(self):
        return model_to_dict(self, fields=[field.name for field in
                             self._meta.fields])
