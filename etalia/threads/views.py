# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic import UpdateView, FormView, DetailView, CreateView
from braces.views import LoginRequiredMixin
from django.core.urlresolvers import reverse

from etalia.core.mixins import AjaxableResponseMixin
from .forms import ThreadCreateForm
from .models import Thread


class ThreadCreate(LoginRequiredMixin, AjaxableResponseMixin, CreateView):

    template_name = 'threads/thread_create.html'
    form_class = ThreadCreateForm

    def get_success_url(self, **kwargs):
        return reverse('threads:thread', kwargs={'pk': self.object.id})

    def get_initial(self):
        return {
            'user_id': self.request.user.id,
            'papers': self.request.user.lib.papers.all()
        }

    def form_valid(self, form):
        form.instance.user = self.request.user
        return super(ThreadCreate, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        return {'redirect': self.get_success_url()}


create = ThreadCreate.as_view()


class ThreadView(LoginRequiredMixin, DetailView):

    model = Thread
    template_name = 'threads/thread.html'

    def get_context_data(self, **kwargs):
        context = super(ThreadView, self).get_context_data(**kwargs)
        context['has_joined'] = self.object.has_joined(self.request.user)
        context['members'] = self.object.members.all()
        context['posts'] = self.object.posts.all().select_related('comments')




thread = ThreadView.as_view()