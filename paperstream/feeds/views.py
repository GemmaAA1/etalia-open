# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import random

from django.http import JsonResponse
from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.views.generic import DeleteView, RedirectView, CreateView
from django.shortcuts import get_object_or_404
from django.db.models import Q
from django.forms.utils import ErrorList
from django.core.exceptions import ValidationError

from braces.views import LoginRequiredMixin

from endless_pagination.views import AjaxListView

from paperstream.core.mixins import ModalMixin, AjaxableResponseMixin
from paperstream.library.models import Paper, Stats
from paperstream.users.models import UserTaste

from .models import UserFeed, UserFeedMatchPaper, UserFeedSeedPaper, \
    DiscoverFeedPaper
from .forms import CreateUserFeedForm
from .tasks import update_feed as async_update_feed


class FeedView(LoginRequiredMixin, ModalMixin, AjaxListView):
    """ClassView for displaying a UserFeed instance"""
    model = UserFeedMatchPaper
    template_name = 'feeds/feed.html'
    page_template = 'feeds/feed_sub_page.html'
    first_page = 10
    per_page = 5
    context_object_name = 'ufmp_list'

    def get_queryset(self):
        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user, is_ticked=True)\
            .values('paper')
        # papers_ticked = []

        query_set = UserFeedMatchPaper.objects\
            .filter(feed=self.userfeed)\
            .exclude(paper__in=papers_ticked)\

        query_set = self.filter_queryset(query_set)

        return query_set

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(paper__title__icontains=q) |
                                   Q(paper__abstract__icontains=q) |
                                   Q(paper__journal__title__icontains=q) |
                                   Q(paper__authors__last_name__icontains=q) |
                                   Q(paper__authors__first_name__icontains=q))\
                .distinct()
        # No q is specified so we return queryset
        return queryset

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context,
                       **super(FeedView, self).get_context_data(**kwargs))
        context['userfeed'] = self.userfeed

        context['first_page'] = self.first_page
        context['per_page'] = self.per_page

        # Get user tastes label
        user_taste = UserTaste.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'is_liked', 'is_ticked')
        # reformat to dict
        user_taste = dict((key, {'liked': v1, 'is_ticked': v2})
                          for key, v1, v2 in user_taste)
        context['user_taste'] = user_taste

        # Get library stats
        context['stats'] = Stats.objects.last()

        if self.request.GET.get("query"):
            context['current_query'] = self.request.GET.get("query")

        return context

    def get_userfeed(self, request, **kwargs):
        feed_name = kwargs.get('name')
        userfeed = get_object_or_404(UserFeed,
                                     name=feed_name,
                                     user=request.user)
        return userfeed

    def get(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        return super(FeedView, self).get(request, *args, **kwargs)

feed_view = FeedView.as_view()


class FeedMainView(LoginRequiredMixin, RedirectView):
    """Redirect to main feed"""

    pattern_name = 'feeds:feed'
    permanent = False

    def get_redirect_url(self, *args, **kwargs):
        # Check if userfeed pk matched db and current logged user
        return super(FeedMainView, self).get_redirect_url(name='main')

feed_main = FeedMainView.as_view()


class CreateFeedView(LoginRequiredMixin, AjaxableResponseMixin, CreateView):
    """ClassView to create a new UserFeed instance"""
    model = UserFeed
    form_class = CreateUserFeedForm
    template_name = 'feeds/feed.html'

    def get_success_url(self):
        return reverse('feeds:modify-feed',
                       kwargs={'name': self.object.name})

    def form_valid(self, form):
        self.object = form.save(commit=False)
        self.object.user = self.request.user
        try:
            self.object.full_clean()
        except ValidationError:
            form._errors['name'] = \
                ErrorList(['You already have a feed name like this'])
            return super(CreateFeedView, self).form_invalid(form)

        return super(CreateFeedView, self).form_valid(form)

    def form_invalid(self, form):
        return super(CreateFeedView, self).form_invalid(form)

    def get_ajax_data(self):
        data = {'redirect':
                reverse('feeds:modify-feed', kwargs={'name': self.object.name})
                }
        return data

create_feed_view = CreateFeedView.as_view()


class ModifyFeedView(LoginRequiredMixin, ModalMixin, ListView):
    """ClassView to modify paper seed of user feed"""
    model = Paper
    template_name = 'feeds/feed_settings.html'
    paginate_by = 50

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(title__icontains=q) |
                                   Q(abstract__icontains=q) |
                                   Q(journal__title__icontains=q) |
                                   Q(authors__last_name__icontains=q) |
                                   Q(authors__first_name__icontains=q))\
                .distinct()
        # No q is specified so we return queryset
        return queryset

    def get_queryset(self):
        # get all user lib paper
        queryset = self.request.user.lib.papers.all()
        # filter through search input
        queryset = self.filter_queryset(queryset)
        return queryset

    def get_context_data(self, **kwargs):
        context = super(ModifyFeedView, self).get_context_data(**kwargs)
        context['userfeed'] = self.userfeed
        context['seed_pks'] = self.userfeed.papers_seed.all()\
            .values_list('pk', flat=True)
        return context

    def get_userfeed(self, request, **kwargs):
        feed_name = kwargs.get('name', None)
        if feed_name:
            userfeed = get_object_or_404(UserFeed,
                                         name=feed_name,
                                         user=request.user)
        else:
            userfeed = get_object_or_404(UserFeed,
                                         user=request.user,
                                         name='main')
        return userfeed

    def get(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        return super(ModifyFeedView, self).get(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        seed_pks = self.userfeed.papers_seed.all().values_list('pk', flat=True)
        if request.POST:
            pk = int(request.POST.get('pk'))
            if pk in seed_pks:      # remove it
                UserFeedSeedPaper.objects\
                    .get(feed=self.userfeed, paper_id=pk)\
                    .delete()
                data = {str(pk): ''}
            else:       # add it
                UserFeedSeedPaper.objects\
                    .create(feed=self.userfeed, paper_id=pk)
                data = {str(pk): 'success'}
            if request.is_ajax():
                return JsonResponse(data)
            else:
                return super(ModifyFeedView, self).get(self, request, *args,
                                                       **kwargs)

modify_feed_view = ModifyFeedView.as_view()


class DeleteFeedView(LoginRequiredMixin, ModalMixin, DeleteView):
    """ClassView to delete a UserFeed instance"""

    model = UserFeed
    success_url = reverse_lazy('feeds:main')

delete_feed_view = DeleteFeedView.as_view()


class UpdateFeedView(LoginRequiredMixin, ModalMixin, RedirectView):

    pattern_name = 'feeds:feed'
    permanent = False

    def get_redirect_url(self, *args, **kwargs):
        # Check if userfeed pk matched db and current logged user
        feed_name = kwargs.get('name', '')
        # update feed async
        async_update_feed.delay(self.request.user.pk, feed_name=feed_name)
        return super(UpdateFeedView, self).get_redirect_url(*args, **kwargs)

update_feed_view = UpdateFeedView.as_view()


class DiscoverView(LoginRequiredMixin, ModalMixin, AjaxListView):
    """ClassView for displaying a UserFeed instance"""
    model = DiscoverFeedPaper
    template_name = 'feeds/discover.html'
    page_template = 'feeds/discover_sub_page.html'
    first_page = 30
    per_page = 20
    context_object_name = 'dfp_list'

    def get_queryset(self):
        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user, is_ticked=True)\
            .values('paper')

        query_set = DiscoverFeedPaper.objects\
            .filter(discover_feed=self.request.user.discover)\
            .exclude(paper__in=papers_ticked)

        # filter from searchbox
        query_set = self.filter_queryset(query_set)

        # let's shuffle the results
        return query_set.order_by("?")

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(paper__title__icontains=q) |
                                   Q(paper__abstract__icontains=q) |
                                   Q(paper__journal__title__icontains=q) |
                                   Q(paper__authors__last_name__icontains=q) |
                                   Q(paper__authors__first_name__icontains=q))\
                .distinct()
        # No q is specified so we return queryset
        return queryset

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context,
                       **super(DiscoverView, self).get_context_data(**kwargs))

        context['first_page'] = self.first_page
        context['per_page'] = self.per_page

        # Get user tastes label
        user_taste = UserTaste.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'is_liked', 'is_ticked')
        # reformat to dict
        user_taste = dict((key, {'liked': v1, 'is_ticked': v2})
                          for key, v1, v2 in user_taste)
        context['user_taste'] = user_taste

        # Get library stats
        context['stats'] = Stats.objects.last()

        if self.request.GET.get("query"):
            context['current_query'] = self.request.GET.get("query")

        return context

discover = DiscoverView.as_view()


@login_required
def ajax_user_feed_message(request, name):
    if request.method == 'GET':
        userfeed = get_object_or_404(UserFeed,
                                     name=name,
                                     user=request.user)
        if userfeed.state == 'IDL':
            data = {'done': True,
                    'url': str(reverse_lazy('feeds:feed',
                                        kwargs={'name': userfeed.name}))}
        else:
            data = {'done': False,
                    'message': userfeed.message}
        return JsonResponse(data)