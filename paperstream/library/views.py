from django.shortcuts import render
from django.views.generic.list import ListView
from django.views.generic.detail import DetailView
from django.conf import settings
from django.shortcuts import get_object_or_404
from django.conf import settings

from .models import Journal, Paper
from .constants import PAPER_TYPE


def library(request):

    context = {'journal_count': Journal.objects.count(),
               'paper_count': Paper.objects.count()}
    return render(request, 'library/library.html', context)


class JournalsListView(ListView):
    model = Journal
    template_name = 'library/journals.html'
    paginate_by = settings.ITEMS_PER_PAGE

journals = JournalsListView.as_view()


class JournalView(ListView):
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


class PaperView(DetailView):

    model = Paper
    template_name = 'library/paper.html'

    def get_context_data(self, **kwargs):
        context = super(PaperView, self).get_context_data(**kwargs)
        paper_ = kwargs['object']
        model = self.request.user.settings.model

        # Get neighbors papers
        neigh = paper_.neighbors.get(lsh__model=model, lsh__time_lapse=-1)
        if neigh.neighbors:
            neigh_pk = neigh.neighbors[:settings.NUMBER_OF_NEIGHBORS]
            neigh_fetched_d = dict(
                [(p.pk, p) for p in Paper.objects.filter(pk__in=neigh_pk)])
            context['neighbors'] = \
                [neigh_fetched_d[key] for key in neigh_pk]
        else:
            context['neighbors'] = []

        context['paper_type'] = dict(PAPER_TYPE)[kwargs['object'].type]
        return context

paper = PaperView.as_view()


class PapersListView(ListView):
    model = Paper
    template_name = 'library/papers.html'
    paginate_by = settings.ITEMS_PER_PAGE

papers = PapersListView.as_view()

