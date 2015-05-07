from django.shortcuts import render
from django.views.generic.list import ListView
from django.conf import settings
from .models import Journal, Paper

# Create your views here.


class JournalsListView(ListView):
    model = Journal
    template_name = 'journals.html'
    paginate_by = settings.ITEMS_PER_PAGE

view_journals = JournalsListView.as_view()


class JournalListView(ListView):
    model = Paper
    template_name = 'journal.html'
    paginate_by = settings.ITEMS_PER_PAGE

    def get_context_data(self, **kwargs):
        context = super(JournalListView, self).get_context_data(**kwargs)
        context['journal'] = Journal.objects.get(id=self.kwargs['id'])
        return context

    def get_queryset(self):
        papers = Paper.objects.filter(journal__id=self.kwargs['id'])
        return papers

view_journal = JournalListView.as_view()

def view_library(request):
    return render(request, 'library.html')

def view_papers(request):
    pass


def view_paper(request, paper_id):
    pass
