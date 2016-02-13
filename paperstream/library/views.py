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

from paperstream.nlp.models import PaperNeighbors, Model, MostSimilar
from paperstream.core.mixins import ModalMixin
from paperstream.users.models import UserTaste
from .models import Journal, Paper
from .constants import PAPER_TYPE


def library(request):

    context = {'journal_count': Journal.objects.count(),
               'paper_count': Paper.objects.count()}
    return render(request, 'library/library.html', context)


class JournalsListView(ModalMixin, ListView):
    model = Journal
    template_name = 'library/journals.html'
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(JournalsListView, self).get_context_data(**kwargs)
        context['journals_count'] = Journal.objects.count()
        return context

journals = JournalsListView.as_view()


class JournalView(ModalMixin, ListView):
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


class PaperView(ModalMixin, DetailView):

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
                ut = UserTaste.objects.get(user=self.request.user, paper=paper_)
                context['is_pinned'] = ut.is_pinned
            except UserTaste.DoesNotExist:
                pass

            if paper_ in self.request.user.lib.papers.all():
                context['is_in_lib'] = True
            else:
                context['is_in_lib'] = False

        return context

paper_slug = PaperView.as_view()


class PaperViewPk(RedirectView):
    """Redirect to slug paper url"""

    permanent = False
    query_string = True

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        return paper.get_absolute_url() + '?' + self.request.GET.urlencode()

paper = PaperViewPk.as_view()


class PaperNeighborsView(LoginRequiredMixin, ListView):

    model = Paper
    template_name = 'library/paper_neighbors.html'
    paper_id = None
    time_span = 30

    def get_queryset(self):

        if self.request.is_ajax():
            data = self.request.GET.dict()
            self.time_span = int(data.get('time_span', self.time_span)) \
                             or self.time_span
            self.paper_id = int(data.get('id'))

            paper_ = Paper.objects.get(id=self.paper_id)

            # Get active MostSimilar
            ms = MostSimilar.objects.filter(is_active=True)

            # Get stored neighbors matches
            try:
                neigh_data = paper_.neighbors.get(ms=ms, time_lapse=self.time_span)
                if not neigh_data.neighbors or not max(neigh_data.neighbors):  # neighbors can be filled with zero if none was found previously
                    neigh_data.delete()
                    raise PaperNeighbors.DoesNotExist
                elif neigh_data.modified > (timezone.now() - timezone.timedelta(days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
                    neighbors = neigh_data.neighbors
                else:
                    raise PaperNeighbors.DoesNotExist
            except PaperNeighbors.DoesNotExist:   # refresh
                try:
                    ms_task = app.tasks['paperstream.nlp.tasks.mostsimilar_{name}'.format(name=self.request.user.settings.stream_model.name)]
                    res = ms_task.apply_async(args=('populate_neighbors',),
                                              kwargs={'paper_pk': self.paper_id,
                                                      'time_lapse': self.time_span},
                                              timeout=10,
                                              soft_timeout=5)
                    neighbors = res.get()
                except SoftTimeLimitExceeded:
                    neighbors = []
                except KeyError:
                    raise

            neigh_pk_list = [neigh for neigh in neighbors[:settings.NUMBER_OF_NEIGHBORS] if neigh]

            clauses = ' '.join(['WHEN id=%s THEN %s' % (pk, i)
                                for i, pk in enumerate(neigh_pk_list)])
            ordering = 'CASE %s END' % clauses
            return Paper.objects.filter(pk__in=neigh_pk_list).extra(
               select={'ordering': ordering}, order_by=('ordering',))

    def get_context_data(self, **kwargs):
        context = super(PaperNeighborsView, self).get_context_data(**kwargs)

        context['data'] = json.dumps({
            'id': self.paper_id,
            'time-span': self.time_span,
        })
        return context

paper_neighbors = PaperNeighborsView.as_view()