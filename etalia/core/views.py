# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import operator
import json
from functools import reduce
from collections import Counter

from django.shortcuts import render, redirect
from django.db.models import Q
from django.utils import timezone
from django.conf import settings

from endless_pagination.views import AjaxListView

from etalia.altmetric.models import AltmetricModel
from etalia.users.models import Author
from etalia.library.models import PaperUser
from etalia.library.constants import PAPER_PINNED, PAPER_BANNED
from etalia.core.utils import AttrDict
from .tasks import failing_task


def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:stream')
    else:
        # Get some trending altmetric matches
        d = timezone.datetime.now().date() - \
            timezone.timedelta(days=settings.LANDING_ACTIVE_PAPERS_TIME_IN_DAYS)
        ten_recent_most_active_paper = AltmetricModel.objects\
            .filter(Q(paper__date_ep__gt=d) |
                    (Q(paper__date_ep=None) & Q(paper__date_pp__gt=d)))\
            [:settings.LANDING_ACTIVE_PAPERS_NUMBER]
        papers = [p.paper for p in ten_recent_most_active_paper]
        context = {'active_papers': papers}
        return render(request, 'landing.html', context=context)


def about(request):
    context = {}
    return render(request, 'pages/about.html', context=context)


def terms(request):
    context = {}
    return render(request, 'pages/terms.html', context=context)


def support(request):
    context = {}
    return render(request, 'pages/support.html', context=context)


def help(request):
    context = {}
    return render(request, 'pages/help.html', context=context)


def test(request):
    return render(request, 'test.html', {})


class BasePaperListView(AjaxListView):
    """Abstract View for displaying a feed type instance"""
    first_page = 10
    per_page = 5
    context_object_name = 'object_list'
    size_max_journal_filter = 40
    size_max_author_filter = 40
    return_filter = True
    original_qs = None
    template_name = None
    model = None
    control_session_name = ''

    def __init__(self, *args, **kwargs):
        super(BasePaperListView, self).__init__(*args, **kwargs)
        self.journals_filter = []
        self.authors_filter = []
        self.time_span = 30
        self.like_flag = False
        self.context_settings = None
        self.search_query = ''
        self.cluster = 0
        self.endless_only = False

    def filter_journals(self, queryset):
        """filter journals"""
        subset = []
        for pk in self.journals_filter:
            subset.append(Q(paper__journal__pk=pk))
        if subset:
            queryset = queryset\
                .filter(reduce(operator.or_, subset))\
                .distinct()

        return queryset

    def filter_authors(self, queryset):
        """filter authors"""
        subset = []
        for pk in self.authors_filter:
            subset.append(Q(paper__authors__pk=pk))
        if subset:
            queryset = queryset\
                .filter(reduce(operator.or_, subset))\
                .distinct()

        return queryset

    def filter_pin(self, queryset):
        like_pks = PaperUser.objects\
            .filter(user=self.request.user,
                    watch=PAPER_PINNED)\
            .values('paper__pk')
        queryset = queryset.filter(paper_id__in=like_pks)

        return queryset

    def filter_search_query(self, queryset):
        """filter search query"""
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

    def filter_time_span(self, queryset):
        """filter queryset based on time span"""
        from_date = (timezone.now() - timezone.timedelta(days=self.time_span)).date()
        return queryset.filter(date__gt=from_date)

    def get_context_endless(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context['first_page'] = self.first_page
        context['per_page'] = self.per_page
        return context

    def get_context_stats(self):
        if self.object_list:
            return {'number_of_papers': self.object_list.count()}
        else:
            return {'number_of_papers': 0}

    def get_context_usertaste(self):
        user_taste = PaperUser.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'watch')
        # reformat to dict
        user_taste = dict((key, {'is_pinned': w == PAPER_PINNED,
                                 'is_banned': w == PAPER_BANNED})
                          for key, w in user_taste)
        return {'user_taste': user_taste}

    def get_context_userlib(self):
        return {'user_lib':
                    self.request.user.lib.papers
                        .all()
                        .values_list('pk', flat=True)}

    def get_context_journals_filter(self):

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

        return {'journals': journals}

    def get_context_authors_filter(self):

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
            authors.append((auth.pk, auth.print_full))

        return {'authors': authors}

    def get_context_search_query(self):
        return {'search_query': self.search_query}

    def get_context_time_span(self):
        return {'time_span': self.time_span}

    def get_context_new_objects_since_last_login(self):
        return {
            'new_papers':
                self.object_list.filter(new=True).values_list('paper_id', flat=True)
        }

    def get_context_control_session(self):
        if self.request.session.get(self.control_session_name):
            return {
                'control_session':
                    AttrDict(self.request.session[self.control_session_name])
            }
        else:
            return {
                'control_session': None
            }

    def get_input_data(self):
        """Get ajax args"""

        # Get data from ajax call
        if self.request.is_ajax():
            if self.request.GET.dict().get('data'):
                data = json.loads(self.request.GET.dict().get('data'))
                #
                self.populate_attr(data)
                # store controls data in session
                self.store_controls_in_session(data)
        # Get data from session
        elif self.request.session.get(self.control_session_name):
            data = self.request.session[self.control_session_name]
            self.populate_attr(data)
        # Get default data
        else:
            data = {
                'time_span': self.time_span,
                'cluster': self.cluster,
                'pin': self.like_flag,
                'search_query': self.search_query,
                'filters': [{'id': 'journal', 'pk': []}, {'id': 'author', 'pk': []}],
            }
            # store controls data in session
            self.store_controls_in_session(data)

    def populate_attr(self, data):
        try:
            self.time_span = int(data.get('time_span', self.time_span) or self.time_span)
            self.cluster = int(data.get('cluster', self.cluster) or self.cluster)
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

    def store_controls_in_session(self, data):
        """Store controls (filter, pin, time-span to session) in session"""
        self.request.session[self.control_session_name] = data

    def get_context_settings(self):
        raise NotImplemented


def test_failing_task(request):
    failing_task.delay()
    return redirect('feeds:stream')




