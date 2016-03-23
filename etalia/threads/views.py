# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic import UpdateView, FormView, DetailView, CreateView

from braces.views import LoginRequiredMixin
from .forms import ThreadCreateForm


class ThreadCreate(LoginRequiredMixin, FormView):

    template_name = 'threads/thread_create.html'
    form_class = ThreadCreateForm

    def get_initial(self):
        return {'author': self.request.user,
                'paper': self.request.user.lib.papers.all()}

    #
    # def get_success_url(self):
    #     pass


create = ThreadCreate.as_view()