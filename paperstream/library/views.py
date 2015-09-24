# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import render
from django.views.generic.list import ListView
from django.views.generic.detail import DetailView
from django.shortcuts import get_object_or_404
from django.conf import settings
from django.utils import timezone

from config.celery import celery_app as app

from paperstream.nlp.models import PaperNeighbors, Model
from paperstream.core.mixins import ModalMixin
from .models import Journal, Paper
from .constants import PAPER_TYPE


def library(request):

    context = {'journal_count': Journal.objects.count(),
               'paper_count': Paper.objects.count()}
    return render(request, 'library/library.html', context)


class JournalsListView(ModalMixin, ListView):
    model = Journal
    template_name = 'library/journals.html'
    paginate_by = settings.ITEMS_PER_PAGE

journals = JournalsListView.as_view()


class JournalView(ModalMixin, ListView):
    # Journal view display a list of papers from the journal
    model = Paper
    template_name = 'library/journal.html'
    paginate_by = settings.ITEMS_PER_PAGE

    def get_context_data(self, **kwargs):
        context = super(JournalView, self).get_context_data(**kwargs)
        context['journal'] = get_object_or_404(Journal, id=self.kwargs['pk'])
        return context

    def get_queryset(self):
        papers_in_jou = Paper.objects.filter(journal__id=self.kwargs['pk'])
        return papers_in_jou

journal = JournalView.as_view()


class PaperView(ModalMixin, DetailView):

    model = Paper
    template_name = 'library/paper.html'

    def get_context_data(self, **kwargs):
        context = super(PaperView, self).get_context_data(**kwargs)
        paper_ = kwargs['object']
        if self.request.user.is_authenticated():
            model = self.request.user.settings.model
        else:
            model = Model.objects.first()

        # Get stored neighbors papers
        try:
            neigh_data = paper_.neighbors.get(lsh__model=model, lsh__time_lapse=-1)
            if neigh_data.modified > (timezone.now() - timezone.timedelta(days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
                neighbors = neigh_data.neighbors
            else:
                raise PaperNeighbors.DoesNotExist
        except PaperNeighbors.DoesNotExist:   # refresh
            try:
                lsh_task = app.tasks['paperstream.nlp.tasks.lsh_{name}_{time_lapse}'.format(
                    name=model.name,
                    time_lapse=-1)]
                res = lsh_task.delay(paper_.pk, 'populate_neighbors')
                neighbors = res.get()
            except KeyError:
                raise

        neigh_pk = neighbors[:settings.NUMBER_OF_NEIGHBORS]
        neigh_fetched_d = dict(
            [(p.pk, p) for p in Paper.objects.filter(pk__in=neigh_pk)])
        context['neighbors'] = \
            [neigh_fetched_d[key] for key in neigh_pk]

        context['paper_type'] = dict(PAPER_TYPE)[kwargs['object'].type]
        return context

paper = PaperView.as_view()


class PapersListView(ModalMixin, ListView):
    model = Paper
    template_name = 'library/papers.html'
    paginate_by = settings.ITEMS_PER_PAGE

papers = PapersListView.as_view()

