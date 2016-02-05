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
from django.views.generic import DeleteView, CreateView
from django.shortcuts import get_object_or_404, redirect
from django.db.models import Q
from django.forms.utils import ErrorList
from django.core.exceptions import ValidationError
from django.conf import settings
from django.utils import timezone


from braces.views import LoginRequiredMixin

from endless_pagination.views import AjaxListView

from paperstream.core.mixins import ModalMixin, AjaxableResponseMixin
from paperstream.library.models import Paper, Stats, Author
from paperstream.users.models import UserTaste

from .models import Stream, StreamMatches, StreamSeeds, \
    TrendMatches
from .forms import CreateUserFeedForm
from .tasks import update_stream, update_trend, reset_stream, reset_trend


class BasePaperListView(LoginRequiredMixin, AjaxListView):
    """Abstract View for displaying a feed type instance"""
    first_page = 10
    per_page = 5
    context_object_name = 'object_list'
    size_max_journal_filter = 40
    size_load_journal_filter = 10
    size_max_author_filter = 40
    size_load_author_filter = 10
    original_qs = None

    def __init__(self, *args, **kwargs):
        super(BasePaperListView, self).__init__(*args, **kwargs)
        self.journals_filter = []
        self.authors_filter = []
        self.time_span = 7
        self.like_flag = False
        self.context_settings = None
        self.search_query = ''
        self.cluster = None

    def get_context_settings(self):
        raise NotImplemented

    def filter_queryset(self, queryset):

        # filter journals
        if self.journals_filter:
            subset = []
            for pk in self.journals_filter:
                subset.append(Q(paper__journal__pk=pk))
            if subset:
                queryset = queryset\
                    .filter(reduce(operator.or_, subset))\
                    .distinct()

        # filter authors
        if self.authors_filter:
            subset = []
            for pk in self.authors_filter:
                subset.append(Q(paper__authors__pk=pk))
            if subset:
                queryset = queryset\
                    .filter(reduce(operator.or_, subset))\
                    .distinct()

        # filter pinned
        if self.like_flag:
            like_pks = UserTaste.objects\
                .filter(user=self.request.user,
                        is_liked=True)\
                .values('paper__pk')
            queryset = queryset.filter(paper_id__in=like_pks)

        # search query
        if self.search_query:
            # return a filtered queryset
            queryset = queryset.filter(Q(paper__title__icontains=self.search_query) |
                                   Q(paper__abstract__icontains=self.search_query) |
                                   Q(paper__journal__title__icontains=self.search_query) |
                                   Q(paper__authors__last_name__icontains=self.search_query) |
                                   Q(paper__authors__first_name__icontains=self.search_query))\
                .distinct()

        return queryset

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context, **super(BasePaperListView, self).get_context_data(**kwargs))

        context['first_page'] = self.first_page
        context['per_page'] = self.per_page
        context['number_of_papers'] = self.object_list.count()

        # Get user tastes label
        user_taste = UserTaste.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'is_liked', 'is_ticked')
        # reformat to dict
        user_taste = dict((key, {'liked': v1, 'is_ticked': v2})
                          for key, v1, v2 in user_taste)
        context['user_taste'] = user_taste

        # Get paper is_in_lib
        context['user_lib'] = self.request.user.lib.papers.all()\
            .values_list('pk', flat=True)

        # BUILD FILTERS
        # Journal filter
        # Retrieve stream journal data for filters
        journals = self.original_qs\
            .values_list('paper__journal__pk',
                         'paper__journal__title')

        # Retrieve journal list from user lib
        journals_userlib = self.request.user.lib.userlibjournal_set.all()\
            .values_list('journal__pk',
                         'journal__title')

        # Count journal occurrence in current stream
        journals_counter = Counter(journals)
        journals_occ_sorted = sorted(journals_counter, key=journals_counter.get,
                                     reverse=True)[:self.size_max_journal_filter]

        # Concatenate userlibjournal and journal in current stream,
        # checking for uniqueness
        journal_filter = []
        for j in journals_userlib:
            if j in journals_occ_sorted:
                journal_filter.append(j)
        for j in journals_occ_sorted:
            if j not in journal_filter:
                journal_filter.append(j)

        # group journals by block of size size_load_journal_filter
        block = 0
        context['journals'] = []
        for i, j in enumerate(journal_filter):
            if not i % self.size_load_journal_filter:
                block += 1
            # (journal.pk, journal.title, occurence, block )
            context['journals'].append((j[0], j[1], journals_counter[j], block))
        if len(context['journals']) <= self.size_load_journal_filter:
            context['journals_show_more'] = False
        else:
            context['journals_show_more'] = True

        # Author filter
        # Retrieve stream authors data for filters
        authors = self.original_qs\
            .values_list('paper__authors', flat=True)\
            .distinct()

        # retrieve user lib authors and count occurence
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
        # build the query ordered based on author occurrence in user library
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
            context['authors'].append((auth.pk, auth.print_full, block))
        if len(context['authors']) <= self.size_load_author_filter:
            context['authors_show_more'] = False
        else:
            context['authors_show_more'] = True

        # Set filter context
        context['search_query'] = self.search_query

        if self.journals_filter:
            context['journals_filter'] = self.journals_filter
        else:
            context['journals_filter'] = []

        if self.authors_filter:
            context['authors_filter'] = [pk for pk in self.authors_filter]
        else:
            context['authors_filter'] = []

        context['time_span'] = self.time_span

        return context

    def update_args(self):
        """Retrieve ajax and GET arguments"""

        get_args = self.request.GET.dict()

        # Get time-span GET arg
        try:
            self.time_span = int(get_args.pop('time-span'))
        except KeyError:
            pass

        # What remains should be ajax args
        if self.request.is_ajax():
            try:
                data = json.loads(list(get_args.keys())[0])
                if data.get('source') == 'filter':
                    # get data
                    self.journals_filter = data.get('journals')
                    self.authors_filter = data.get('authors')
                    self.like_flag = data.get('like_flag')
                    self.search_query = data.get('search_query')
                    self.time_span = data.get('time_span') or self.time_span
            except ValueError:  # likely data from AjaxListView
                pass

    def trim_time_span(self, queryset):
        """Trim queryset based on time span"""
        from_date = (timezone.now() - timezone.timedelta(days=self.time_span)).date()
        return queryset.filter(date__gt=from_date)


class StreamView(BasePaperListView):
    # TODO: should be depreciated to StreamView2
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

    def get_queryset(self):

        # get ticked/rejected paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(stream__name=self.kwargs.get('name', 'main'),
                    stream__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

        # Retrieve get args
        self.update_args()

        # trim time span
        self.original_qs = self.trim_time_span(self.original_qs)

        # filter
        query_set = self.filter_queryset(self.original_qs)

        return query_set

stream_view = StreamView.as_view()


class TrendView(BasePaperListView):
    # TODO: should be depreciated for TrendView2
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

    def get_queryset(self):

        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(trend__name=self.kwargs.get('name', 'main'),
                    trend__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

        # Retrieve get args
        self.update_args()

        # trim time span
        self.original_qs = self.trim_time_span(self.original_qs)

        # filter
        query_set = self.filter_queryset(self.original_qs)

        return query_set

trend_view = TrendView.as_view()


class BasePaperListView2(LoginRequiredMixin, AjaxListView):
    """Abstract View for displaying a feed type instance"""
    first_page = 10
    per_page = 5
    context_object_name = 'object_list'
    size_max_journal_filter = 40
    size_load_journal_filter = 10
    size_max_author_filter = 40
    size_load_author_filter = 10
    original_qs = None
    return_filter = True

    def __init__(self, *args, **kwargs):
        super(BasePaperListView2, self).__init__(*args, **kwargs)
        self.journals_filter = []
        self.authors_filter = []
        self.time_span = 30
        self.like_flag = False
        self.context_settings = None
        self.search_query = ''
        self.cluster = None
        self.endless_only = False

    def get_context_settings(self):
        raise NotImplemented

    def filter_queryset(self, queryset):

        # filter journals
        if self.journals_filter:
            subset = []
            for pk in self.journals_filter:
                subset.append(Q(paper__journal__pk=pk))
            if subset:
                queryset = queryset\
                    .filter(reduce(operator.or_, subset))\
                    .distinct()

        # filter authors
        if self.authors_filter:
            subset = []
            for pk in self.authors_filter:
                subset.append(Q(paper__authors__pk=pk))
            if subset:
                queryset = queryset\
                    .filter(reduce(operator.or_, subset))\
                    .distinct()

        # filter pinned
        if self.like_flag:
            like_pks = UserTaste.objects\
                .filter(user=self.request.user,
                        is_liked=True)\
                .values('paper__pk')
            queryset = queryset.filter(paper_id__in=like_pks)

        # search query
        if self.search_query:
            subset = []
            for word in self.search_query.split():
                subset.append(Q(paper__title__icontains=word) |
                              Q(paper__abstract__icontains=word) |
                              Q(paper__journal__title__icontains=word) |
                              Q(paper__authors__last_name__icontains=word) |
                              Q(paper__authors__first_name__icontains=word))
            if subset:
                queryset = queryset\
                    .filter(reduce(operator.and_, subset))\
                    .distinct()

        return queryset

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context, **super(BasePaperListView2, self).get_context_data(**kwargs))

        context['first_page'] = self.first_page
        context['per_page'] = self.per_page
        context['number_of_papers'] = self.object_list.count()

        # Get user tastes label
        user_taste = UserTaste.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'is_liked', 'is_ticked')
        # reformat to dict
        user_taste = dict((key, {'liked': v1, 'is_ticked': v2})
                          for key, v1, v2 in user_taste)
        context['user_taste'] = user_taste

        # Get paper is_in_lib
        context['user_lib'] = self.request.user.lib.papers.all()\
            .values_list('pk', flat=True)

        if self.return_filter:
            context = dict(context, **self.get_context_filter())

        return context

    def get_context_filter(self):
        """Build context filter"""
        context = {}

        # Journal filter
        # Retrieve stream journal data for filters
        journals = self.original_qs\
            .values_list('paper__journal__pk',
                         'paper__journal__title')

        # Retrieve journal list from user lib
        journals_userlib = self.request.user.lib.userlibjournal_set.all()\
            .values_list('journal__pk',
                         'journal__title')

        # Count journal occurrence in current stream
        journals_counter = Counter(journals)
        journals_occ_sorted = sorted(journals_counter, key=journals_counter.get,
                                     reverse=True)[:self.size_max_journal_filter]

        # Concatenate userlibjournal and journal in current stream,
        # checking for uniqueness
        journal_filter = []
        for j in journals_userlib:
            if j in journals_occ_sorted:
                journal_filter.append(j)
        for j in journals_occ_sorted:
            if j not in journal_filter:
                journal_filter.append(j)

        # group journals by block of size size_load_journal_filter
        journals = []
        for i, j in enumerate(journal_filter):
            journals.append((j[0], j[1], journals_counter[j]))

        # Author filter
        # Retrieve stream authors data for filters
        authors = self.original_qs\
            .values_list('paper__authors', flat=True)\
            .distinct()

        # retrieve user lib authors and count occurence
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
        # build the query ordered based on author occurrence in user library
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
        # block = 0
        authors = []
        for i, auth in enumerate(authors_ordered):
            authors.append((auth.pk, auth.print_full, None))

        # Set filter context
        context['search_query'] = self.search_query

        context['time_span'] = self.time_span

        # build JSON object
        filter_ = {
            'filters': [],
            'pin': self.like_flag,
            'time_span': self.time_span,
            'cluster': None,
            'search_query': self.search_query,
        }

        # add to journal filter context if journal is checked
        entries = []
        for j in journals:
            is_checked = j[0] in self.journals_filter
            entries.append({
                'pk': j[0],
                'name': j[1],
                'count': j[2],
                'is_checked': is_checked
            })
        filter_['filters'].append({'id': 'journal', 'entries': entries})

        # add to author filter context if author is checked
        entries = []
        for a in authors:
            is_checked = a[0] in self.authors_filter
            entries.append({
                'pk': a[0],
                'name': a[1],
                'count': None,
                'is_checked': is_checked
            })
        filter_['filters'].append({'id': 'author', 'entries': entries})

        context['filter'] = json.dumps(filter_)

        context['authors'] = []
        block = 1
        for i, auth in enumerate(authors_ordered):
            if not i % 10:
                block += 1
            context['authors'].append((auth.pk, auth.print_full, block))

        context['journals'] = []
        block = 1
        for i, j in enumerate(journal_filter):
            if not i % self.size_load_journal_filter:
                block += 1
            # (journal.pk, journal.title, occurence, block )
            context['journals'].append((j[0], j[1], journals_counter[j], block))

        return context

    def update_args(self):
        """Get ajax args"""

        if self.request.is_ajax():
            if self.request.GET.dict().get('data'):
                data = json.loads(self.request.GET.dict().get('data'))
                try:
                    self.time_span = int(data.get('time_span', self.time_span))
                    self.cluster = int(data.get('cluster', self.cluster))
                    self.like_flag = data.get('pin', self.like_flag)
                    self.search_query = data.get('search_query', '')
                    filters_ = data.get('filters', [])
                    for filter_ in filters_:
                        if filter_.get('id') == 'journal':
                            self.journals_filter = filter_.get('pk')
                        if filter_.get('id') == 'author':
                            self.authors_filter = filter_.get('pk')

                except ValueError:
                    pass

    def trim_time_span(self, queryset):
        """Trim queryset based on time span"""
        from_date = (timezone.now() - timezone.timedelta(days=self.time_span)).date()
        return queryset.filter(date__gt=from_date)


class BaseStreamView2(BasePaperListView2):
    model = StreamMatches
    template_name = 'feeds/feed.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.stream_time_lapse,
            'method': self.request.user.settings.stream_method,
            'model': self.request.user.settings.stream_model,
        }
        return self.context_settings

    def get_queryset(self):

        # get ticked/rejected paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(stream__name=self.kwargs.get('name', 'main'),
                    stream__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

        # Retrieve get args
        self.update_args()

        # trim time span
        self.original_qs = self.trim_time_span(self.original_qs)

        # filter
        query_set = self.filter_queryset(self.original_qs)

        return query_set


class StreamView2(BaseStreamView2):
    page_template = 'feeds/feed_sub_page.html'

    def update_args(self):
        if self.request.GET.dict().get('querystring_key'):  # endless scrolling
            self.return_filter = False
        super(StreamView2, self).update_args()

stream_view2 = StreamView2.as_view()


class StreamView2Filter(BaseStreamView2):
    page_template = 'feeds/feed_sub_page2.html'

stream_view2_filter = StreamView2Filter.as_view()


class BaseTrendView2(BasePaperListView2):
    model = TrendMatches
    template_name = 'feeds/trends.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.trend_time_lapse,
            'method': self.request.user.settings.trend_method,
            'model': self.request.user.settings.trend_model,
        }
        return self.context_settings

    def get_queryset(self):

        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(trend__name=self.kwargs.get('name', 'main'),
                    trend__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

        # Retrieve get args
        self.update_args()

        # trim time span
        self.original_qs = self.trim_time_span(self.original_qs)

        # filter
        query_set = self.filter_queryset(self.original_qs)

        return query_set


class TrendView2(BasePaperListView2):
    page_template = 'feeds/trends_sub_page.html'

trend_view2 = TrendView2.as_view()


class TrendView2Filter(BasePaperListView2):
    page_template = 'feeds/trends_sub_page2.html'

trend_view2_filter = TrendView2Filter.as_view()


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
    template_name = 'feeds/../templates/old/feed_settings.html'
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
def reset_stream_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        stream_name = name
        reset_stream.delay(request.user.pk, stream_name=stream_name,
                           restrict_journal=False)
        data = {'display_update_modal': True,
                'message': 'Stream reset launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_stream_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        stream_name = name
        update_stream.delay(request.user.pk, stream_name=stream_name)
        data = {'display_update_modal': True,
                'message': 'Stream update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_trend_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        trend_name = name
        update_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def reset_trend_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        trend_name = name
        reset_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
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