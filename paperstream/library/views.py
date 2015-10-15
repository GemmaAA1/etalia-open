# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from django.shortcuts import render
from django.views.generic import RedirectView
from django.views.generic.list import ListView
from django.views.generic.detail import DetailView
from django.shortcuts import get_object_or_404
from django.conf import settings
from django.utils import timezone
from django.utils.text import slugify
from django.core.urlresolvers import reverse

from config.celery import celery_app as app

from paperstream.nlp.models import PaperNeighbors, Model
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
    time_lapse_map = {'year': 365,
                      'month': 30,
                      'week': 7}

    def get_context_data(self, **kwargs):
        context = super(PaperView, self).get_context_data(**kwargs)
        paper_ = kwargs['object']
        time_lapse = self.time_lapse_map[self.kwargs.get('time_lapse', 'year')]
        if self.request.user.is_authenticated():
            model = self.request.user.settings.model
        else:
            model = Model.objects.first()

        # Get stored neighbors papers
        try:
            neigh_data = paper_.neighbors.get(model=model, time_lapse=time_lapse)
            if not neigh_data.neighbors:
                neigh_data.delete()
                raise PaperNeighbors.DoesNotExist
            elif neigh_data.modified > (timezone.now() - timezone.timedelta(days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
                neighbors = neigh_data.neighbors
            else:
                raise PaperNeighbors.DoesNotExist
        except PaperNeighbors.DoesNotExist:   # refresh
            try:
                ms_task = app.tasks['paperstream.nlp.tasks.mostsimilar_{name}'.format(
                    name=model.name)]
                res = ms_task.apply_async(args=('populate_neighbors',),
                                          kwargs={'paper_pk': paper_.pk,
                                                  'time_lapse': time_lapse},
                                          timeout=3)
                neighbors = res.get()
            except KeyError:
                raise

        neigh_pk = neighbors[:settings.NUMBER_OF_NEIGHBORS]
        neigh_fetched_d = dict(
            [(p.pk, p) for p in Paper.objects.filter(pk__in=neigh_pk)])
        context['neighbors'] = \
            [neigh_fetched_d[key] for key in neigh_pk]

        context['paper_type'] = dict(PAPER_TYPE)[kwargs['object'].type]
        context['time_lapse'] = self.kwargs.get('time_lapse', 'year')

        if not self.request.user.is_anonymous():
            try:
                ut = UserTaste.objects.get(user=self.request.user, paper=paper_)
                context['is_liked'] = ut.is_liked
            except UserTaste.DoesNotExist:
                pass

        if paper_ in self.request.user.lib.papers.all():
            context['is_in_lib'] = True
        else:
            context['is_in_lib'] = False

        return context

paper_slug = PaperView.as_view()


class PaperViewPk(RedirectView):
    """Redirect to slug paper url for SEO"""

    permanent = True

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        return paper.get_absolute_url()

paper = PaperViewPk.as_view()


class PaperViewPkTime(RedirectView):

    permanent = True

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        return reverse('library:paper-slug-time',
                       kwargs={'pk': paper.pk,
                               'slug': slugify(paper.title),
                               'time_lapse': kwargs.get('time_lapse')})

paper_time = PaperViewPkTime.as_view()


class PapersListView(ModalMixin, ListView):
    model = Paper
    template_name = 'library/papers.html'
    paginate_by = settings.ITEMS_PER_PAGE

papers = PapersListView.as_view()