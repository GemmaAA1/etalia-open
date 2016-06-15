# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json

from django.shortcuts import render
from django.views.generic import RedirectView
from django.views.generic.list import ListView
from django.views.generic.detail import DetailView
from django.shortcuts import get_object_or_404
from django.conf import settings
from django.utils import timezone
from django.utils.text import slugify
from django.core.urlresolvers import reverse
from django.db.models.expressions import RawSQL

from config.celery import celery_app as app
from celery.exceptions import SoftTimeLimitExceeded
from braces.views import LoginRequiredMixin

from etalia.nlp.models import PaperNeighbors, Model, PaperEngine
from etalia.users.mixins import ProfileModalFormsMixin
from etalia.feeds.mixins import CreateFeedModalMixin
from etalia.users.models import UserLibPaper
from etalia.library.models import PaperUser
from etalia.library.constants import PAPER_PINNED, PAPER_BANNED
from .models import Journal, Paper
from .constants import PAPER_TYPE


from django.template.response import TemplateResponse


def my_papers(request):
    return TemplateResponse(
        request,
        'papers/list.html',
        {'control_states': json.dumps(request.session.get('library-control-states', {}))}
    )


def library(request):

    context = {'journal_count': Journal.objects.count(),
               'paper_count': Paper.objects.count()}
    return render(request, 'library/library.html', context)


class JournalsListView(ProfileModalFormsMixin, CreateFeedModalMixin, ListView):
    model = Journal
    template_name = 'library/journals.html'
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(JournalsListView, self).get_context_data(**kwargs)
        context['journals_count'] = Journal.objects.count()
        return context

journals = JournalsListView.as_view()


class JournalView(ProfileModalFormsMixin, CreateFeedModalMixin, ListView):
    # Journal view display a list of matches from the journal
    model = Paper
    template_name = 'library/journal.html'
    paginate_by = settings.LIBRARY_ITEMS_PER_PAGE

    def get_context_data(self, **kwargs):
        context = super(JournalView, self).get_context_data(**kwargs)
        context['journal'] = get_object_or_404(Journal,
                                               id=self.kwargs['pk'])
        return context

    def get_queryset(self):
        papers_in_jou = Paper.objects\
            .filter(journal__id=self.kwargs['pk'],
                    is_trusted=True,
                    type='JOU')\
            .annotate(date=RawSQL("LEAST(date_ep, date_fs, date_pp) ", []))\
            .order_by('-date')

        return papers_in_jou

journal_slug = JournalView.as_view()


class JournalViewPk(RedirectView):
    """Redirect to slug journal url"""

    permanent = True

    def get_redirect_url(self, *args, **kwargs):
        journal = Journal.objects.get(pk=kwargs['pk'])
        return journal.get_absolute_url()

journal = JournalViewPk.as_view()


class PaperView(ProfileModalFormsMixin, CreateFeedModalMixin, DetailView):

    model = Paper
    template_name = 'library/paper.html'
    time_lapse_map = {'year': 365,
                      'month': 30,
                      'week': 7}

    def get_template_names(self):
        if not self.request.is_ajax():
            self.template_name = 'library/paper-page.html'

        return super(PaperView, self).get_template_names()

    def get_context_data(self, **kwargs):
        context = super(PaperView, self).get_context_data(**kwargs)
        paper_ = kwargs['object']
        context['paper_type'] = dict(PAPER_TYPE)[kwargs['object'].type]
        context['time_lapse'] = self.request.GET.get('time-span', 'year')

        if not self.request.user.is_anonymous():
            try:
                ut = PaperUser.objects.get(user=self.request.user, paper=paper_)
                context['is_pinned'] = ut.watch == PAPER_PINNED
            except PaperUser.DoesNotExist:
                pass

            try:
                ulp = self.request.user.lib.userlib_paper.get(paper=paper_)
                context['is_in_lib'] = True
                if ulp.is_trashed:
                    context['is_in_trash'] = True
                else:
                    context['is_in_trash'] = False
            except UserLibPaper.DoesNotExist:
                context['is_in_lib'] = False
                context['is_in_trash'] = False
                pass

        return context

paper_slug = PaperView.as_view()


class PaperViewPk(RedirectView):
    """Redirect to slug paper url"""

    permanent = False
    query_string = True
    pattern_name = 'library:paper-slug'

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(paper.title)
        return super(PaperViewPk, self).get_redirect_url(*args, **kwargs)

paper = PaperViewPk.as_view()


class PaperNeighborsView(LoginRequiredMixin, ListView):

    model = Paper
    template_name = 'library/paper_neighbors.html'
    paper_id = None
    time_span = 30

    def get_queryset(self, **kwargs):
        from .tasks import get_neighbors_papers
        if self.request.is_ajax():
            if self.request.GET.dict().get('data'):
                data = json.loads(self.request.GET.dict().get('data'))
                self.time_span = int(data.get('time_span', self.time_span)) \
                                 or self.time_span

                self.paper_id = int(self.kwargs['pk'])

                return get_neighbors_papers(self.paper_id, self.time_span)

    def get_context_usertaste(self):
        user_taste = PaperUser.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'watch')
        # reformat to dict
        user_taste = dict((key, {'liked': w == PAPER_PINNED,
                                 'is_banned': w == PAPER_BANNED})
                          for key, w in user_taste)
        return {'user_taste': user_taste}

    def get_context_data(self, **kwargs):
        context = super(PaperNeighborsView, self).get_context_data(**kwargs)
        context.update(self.get_context_usertaste())
        context['data'] = json.dumps({
            'time_span': self.time_span,
        })
        return context

paper_neighbors = PaperNeighborsView.as_view()
