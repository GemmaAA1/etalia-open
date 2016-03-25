# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic.base import ContextMixin

from .forms import ThreadPostForm, ThreadPostComment


class ThreadFormsMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(ThreadFormsMixin, self).get_context_data(**kwargs)
        if not self.request.user.is_anonymous():
            context['form_post'] = ThreadPostForm()
            context['form_comment'] = ThreadPostComment()
        return context
