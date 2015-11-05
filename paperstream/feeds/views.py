# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import operator
import json
from functools import reduce
from collections import Counter

from django.http import JsonResponse
from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.views.generic import DeleteView, RedirectView, CreateView
from django.shortcuts import get_object_or_404, redirect
from django.db.models import Q
from django.forms.utils import ErrorList
from django.core.exceptions import ValidationError
from django.conf import settings
from django.db.models.expressions import RawSQL


from braces.views import LoginRequiredMixin

from endless_pagination.views import AjaxListView

from paperstream.core.mixins import ModalMixin, AjaxableResponseMixin
from paperstream.library.models import Paper, Stats, Author
from paperstream.users.models import UserTaste, FeedLayout

from .models import Stream, StreamMatches, StreamSeeds, \
    TrendMatches
from .forms import CreateUserFeedForm
from .tasks import update_stream, update_trend


class BaseFeedView(LoginRequiredMixin, ModalMixin, AjaxListView):
    """Abstract View for displaying a feed type instance"""
    first_page = 10
    per_page = 5
    context_object_name = 'object_list'
    size_max_journal_filter = 40
    size_load_journal_filter = 10
    size_max_author_filter = 40
    size_load_author_filter = 10
    original_qs = None
    journals_filter = None
    journals_filter_flag = 'all'
    authors_filter = None
    authors_filter_flag = 'all'
    sorting_flag = 'relevant'
    like_flag = False
    context_settings = None

    def update_from_filter(self):
        raise NotImplemented

    def get_context_settings(self):
        raise NotImplemented

    def filter_queryset(self, queryset):

        # load stored filter
        self.update_from_filter()

        # journals
        if self.journals_filter_flag == 'all':
            if self.journals_filter:
                subset = []
                for journal, val in self.journals_filter:
                    if not val:
                        subset.append(Q(paper__journal__title__icontains=journal))
                if subset:
                    queryset = queryset.exclude(reduce(operator.or_, subset)).distinct()
        if self.journals_filter_flag == 'none':
            if self.journals_filter:
                subset = []
                for journal, val in self.journals_filter:
                    if val:
                        subset.append(Q(paper__journal__title__icontains=journal))
                if subset:
                    queryset = queryset.filter(reduce(operator.or_, subset)).distinct()
                else:
                    queryset = None
            else:
                queryset = None

        # authors
        if self.authors_filter_flag == 'all':
            if self.authors_filter:
                subset = []
                for pk, val in self.authors_filter:
                    if not val:
                        subset.append(Q(paper__authors__pk=pk))
                if subset:
                    queryset = queryset.exclude(reduce(operator.or_, subset)).distinct()
        if self.authors_filter_flag == 'none':
            if self.authors_filter:
                subset = []
                for pk, val in self.authors_filter:
                    if val:
                        subset.append(Q(paper__authors__pk=pk))
                if subset:
                    queryset = queryset.filter(reduce(operator.or_, subset)).distinct()
                else:
                    queryset = None
            else:
                queryset = None

        # like filter
        if self.like_flag:
            like_pks = UserTaste.objects\
                .filter(user=self.request.user,
                        is_liked=True)\
                .values('paper__pk')
            queryset = queryset.filter(paper_id__in=like_pks)

        # search query
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            queryset = queryset.filter(Q(paper__title__icontains=q) |
                                   Q(paper__abstract__icontains=q) |
                                   Q(paper__journal__title__icontains=q) |
                                   Q(paper__authors__last_name__icontains=q) |
                                   Q(paper__authors__first_name__icontains=q))\
                .distinct()

        # sort queryset
        if self.sorting_flag == 'relevant':
            # this is by default for feed
            pass
        elif self.sorting_flag == 'recent':
            # order by date
            queryset = queryset\
                .annotate(date=RawSQL("SELECT LEAST(date_ep, date_fs, date_pp) "
                                      "FROM library_paper "
                                      "WHERE id = paper_id", []))
            queryset = queryset.order_by('-date')
        elif self.sorting_flag == 'trendy':
            # order by altmetric score
            queryset = queryset.order_by('-paper__altmetric__score')
        elif self.sorting_flag == 'nothing':
            # let's shuffle the results
            queryset = queryset.order_by("?")

        return queryset

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context, **super(BaseFeedView, self).get_context_data(**kwargs))

        user_settings = self.get_context_settings()

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

        # Get paper is_in_lib
        context['user_lib'] = self.request.user.lib.papers.all().values_list('pk', flat=True)

        # Get library stats
        context['stats'] = Stats.objects.last()

        # Get data for filter
        data = self.original_qs.values('paper',
                                       'paper__journal__title',
                                       'paper__authors')
        j_titles = []
        authors = []
        check_papers = []
        for d in data:
            if d['paper'] not in check_papers:  # rows are for different authors
                j_titles.append(d['paper__journal__title'])
                check_papers.append(d['paper'])
            authors.append(d['paper__authors'])

        # journal filter
        journals_counter = Counter(j_titles)
        j_titles = sorted(journals_counter, key=journals_counter.get,
                          reverse=True)[:self.size_max_journal_filter]
        # group authors by block identified with block tag
        block = 0
        context['journals'] = []
        for i, title in enumerate(j_titles):
            if not i % self.size_load_journal_filter:
                block += 1
            context['journals'].append((title, journals_counter[title], block))
        if len(context['journals']) <= self.size_load_journal_filter:
            context['journals_show_more'] = False
        else:
            context['journals_show_more'] = True

        # author filter
        # get user lib authors
        authors_user_lib = self.request.user.lib.papers.all()\
            .values_list('authors', flat=True)
        authors_user_lib_counter = Counter(authors_user_lib)
        authors_user_lib_pk_sorted = sorted(authors_user_lib_counter,
                                            key=authors_user_lib_counter.get,
                                            reverse=True)
        authors_dict = {}
        for auth in authors:
            if auth in authors_user_lib_pk_sorted:
                authors_dict[auth] = authors_user_lib_counter[auth]
            else:
                authors_dict[auth] = 0
        authors_pks = sorted(authors_dict, key=authors_dict.get,
                                 reverse=True)
        clauses = ' '.join(['WHEN id=%s THEN %s' %
                            (pk, i) for i, pk in enumerate(authors_pks)])
        ordering = 'CASE %s END' % clauses
        authors_ordered = Author.objects\
            .filter(pk__in=authors_pks)\
            .exclude(Q(first_name='') & Q(last_name=''))\
            .extra(select={'ordering': ordering},
                   order_by=('ordering',))[:self.size_max_author_filter]
        # group authors by block identified with block tag
        block = 0
        context['authors'] = []
        for i, auth in enumerate(authors_ordered):
            if not i % 10:
                block += 1
            context['authors'].append((auth, auth.print_full, block))
        if len(context['authors']) <= self.size_load_author_filter:
            context['authors_show_more'] = False
        else:
            context['authors_show_more'] = True

        # Set filter context
        if self.request.GET.get('query'):
            context['get_query'] = self.request.GET.get('query')

        if self.journals_filter:
            context['journals_filter'] = [title.lower() for title, val in self.journals_filter if val]
        else:
            context['journals_filter'] = [title.lower() for title in j_titles]
        if self.authors_filter:
            context['authors_filter'] = [pk for pk, val in self.authors_filter if val]
        else:
            context['authors_filter'] = authors_pks

        context['sorting_flag'] = self.sorting_flag

        # time lapse settings
        context['time_lapse'] = \
            dict(settings.NLP_TIME_LAPSE_CHOICES)\
                .get(user_settings.get('time_lapse'))

        return context


class StreamView(BaseFeedView):
    """ClassView for displaying a UserFeed instance"""
    model = StreamMatches
    template_name = 'feeds/feed.html'
    page_template = 'feeds/feed_sub_page.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.stream_time_lapse,
            'method': self.request.user.settings.stream_method,
            'model': self.request.user.settings.stream_model,
        }
        return self.context_settings


    def update_from_filter(self):
        ufl, new = FeedLayout.objects.get_or_create(user=self.request.user)
        if new:
            ufl.stream_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.trend_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.library_filter = {'journals_flag': self.journals_filter_flag,
                                  'authors_flag': self.authors_filter_flag,
                                  'sorting_flag': self.sorting_flag}
            ufl.save()
        else:
            self.journals_filter = ufl.stream_filter.get('journals')
            self.authors_filter = ufl.stream_filter.get('authors')
            self.authors_filter_flag = ufl.stream_filter.get('authors_flag') or self.authors_filter_flag
            self.journals_filter_flag = ufl.stream_filter.get('journals_flag') or self.journals_filter_flag
            self.sorting_flag = ufl.stream_filter.get('sorting_flag') or self.sorting_flag

        # From ajaxable filter
        if self.request.is_ajax():
            try:
                data = json.loads(list(self.request.GET)[0])
                if data['action'] == 'filter':
                    # get data
                    self.journals_filter = data.get('journals')
                    self.authors_filter = data.get('authors')
                    self.authors_filter_flag = data.get('authors_flag') or self.authors_filter_flag
                    self.journals_filter_flag = data.get('journals_flag') or self.journals_filter_flag
                    self.sorting_flag = data.get('sorting_flag') or self.sorting_flag
                    self.like_flag = data.get('like_flag')
                    ufl.stream_filter = {
                        'journals': self.journals_filter,
                        'authors': self.authors_filter,
                        'journals_flag': self.journals_filter_flag,
                        'authors_flag': self.authors_filter_flag,
                        'sorting_flag': self.sorting_flag,
                    }
                    ufl.save()
            except ValueError:  # likely data from AjaxListView
                pass

    def get_queryset(self):

        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(stream__name=self.kwargs.get('name', 'main'),
                    stream__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper', 'paper__journal')

        query_set = self.filter_queryset(self.original_qs)

        return query_set

stream_view = StreamView.as_view()


class TrendView(BaseFeedView):
    """ClassView for displaying a UserFeed instance"""
    model = TrendMatches
    template_name = 'feeds/trends.html'
    page_template = 'feeds/trends_sub_page.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.trend_time_lapse,
            'method': self.request.user.settings.trend_method,
            'model': self.request.user.settings.trend_model,
        }
        return self.context_settings

    def update_from_filter(self):
        ufl, new = FeedLayout.objects.get_or_create(user=self.request.user)
        if new:
            ufl.stream_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.trend_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.library_filter = {'journals_flag': self.journals_filter_flag,
                                  'authors_flag': self.authors_filter_flag,
                                  'sorting_flag': self.sorting_flag}
            ufl.save()
        else:
            self.journals_filter = ufl.trend_filter.get('journals')
            self.authors_filter = ufl.trend_filter.get('authors')
            self.authors_filter_flag = ufl.trend_filter.get('authors_flag') or self.authors_filter_flag
            self.journals_filter_flag = ufl.trend_filter.get('journals_flag') or self.journals_filter_flag
            self.sorting_flag = ufl.trend_filter.get('sorting_flag') or self.sorting_flag

        # From ajaxable filter
        if self.request.is_ajax():
            try:
                data = json.loads(list(self.request.GET)[0])
                if data['action'] == 'filter':
                    # get data
                    self.journals_filter = data.get('journals')
                    self.authors_filter = data.get('authors')
                    self.authors_filter_flag = data.get('authors_flag') or self.authors_filter_flag
                    self.journals_filter_flag = data.get('journals_flag') or self.journals_filter_flag
                    self.sorting_flag = data.get('sorting_flag') or self.sorting_flag
                    self.like_flag = data.get('like_flag')
                    ufl.trend_filter = {
                        'journals': self.journals_filter,
                        'authors': self.authors_filter,
                        'journals_flag': self.journals_filter_flag,
                        'authors_flag': self.authors_filter_flag,
                        'sorting_flag': self.sorting_flag
                    }
                    ufl.save()
            except ValueError:  # likely data from AjaxListView
                pass

    def get_queryset(self):

        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(trend__name=self.kwargs.get('name', 'main'),
                    trend__user=self.request.user)\
            .exclude(paper__in=papers_ticked)

        # filter from searchbox
        query_set = self.filter_queryset(self.original_qs)

        return query_set

trend_view = TrendView.as_view()


class CreateFeedView(LoginRequiredMixin, AjaxableResponseMixin, CreateView):
    """ClassView to create a new UserFeed instance"""
    model = Stream
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
            userfeed = get_object_or_404(Stream,
                                         name=feed_name,
                                         user=request.user)
        else:
            userfeed = get_object_or_404(Stream,
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
                StreamSeeds.objects\
                    .get(feed=self.userfeed, paper_id=pk)\
                    .delete()
                data = {str(pk): ''}
            else:       # add it
                StreamSeeds.objects\
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

    model = Stream
    success_url = reverse_lazy('feeds:main')

delete_feed_view = DeleteFeedView.as_view()


@login_required
def update_stream_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        stream_name = name
        update_stream.delay(request.user.pk, stream_name=stream_name)
        data = {'message': 'Stream update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_trend_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        trend_name = name
        update_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def ajax_user_feed_message(request, name):
    if request.method == 'GET':
        userfeed = get_object_or_404(Stream,
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