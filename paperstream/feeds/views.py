from django.shortcuts import render
from braces.views import LoginRequiredMixin
from django.views.generic.list import ListView

from library.models import Paper

class home_feed(LoginRequiredMixin, ListView):
    model = Paper
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        return Paper.objects.all()

home = home_feed.as_view()