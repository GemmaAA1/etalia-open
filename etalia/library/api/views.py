# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from collections import Counter

from django.views.decorators.cache import never_cache, cache_page
from rest_framework import viewsets, permissions, status, mixins, filters
from rest_framework.decorators import detail_route, list_route
from rest_framework.response import Response
from django.db.models import Q
from django.conf import settings
from django.utils.decorators import method_decorator
from django.shortcuts import get_object_or_404

from etalia.core.api.permissions import IsReadOnlyRequest
from etalia.core.api.mixins import MultiSerializerMixin
from etalia.core.api.filters import DisabledHTMLContextualFilterBackend, \
    DisabledHTMLSearchFilterBackend
from etalia.core.api.permissions import IsOwner, IsSessionAuthenticatedOrReadOnly
from etalia.threads.api.serializers import ThreadSerializer
from etalia.users.constants import USER_INDIVIDUAL
from ..models import Paper, Author, Journal, PaperUser
from ..constants import PAPER_TRASHED
from .serializers import PaperSerializer, JournalSerializer, AuthorSerializer, \
    PaperNestedSerializer, PaperFilterSerializer, PaperUserSerializer, \
    JsonLdSerializer
from ..filters import PaperFilter, MyPaperFilter


class PaperViewSet(MultiSerializerMixin,
                   viewsets.ReadOnlyModelViewSet):
    """Papers

    ### Routes ###

    * **[GET] /papers/**: List of papers
    * **[GET] /papers/<id>**: Detail of paper
    * **[GET] /papers/<id\>/neighbors**: Similar papers
    * **[GET] /papers/<id\>/related-threads**: Related threads
    * **[GET] /papers/<id\>/json-ld**: JSON-LD based on schema.org ScholarlyArticle item

    ### Additional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **doi=(str)**: Fetch by DOI
    * **pmi=(str)**: Fetch by Pub Med Id
    * **arxiv=(str)**: Fetch by Arxiv id
    * **pii=(str)**: Fetch by Elsevier id
    * **title=(str)**: Fetch by title
    * **journal=(str)**: Fetch by journal title
    * **issn=(str)**: Fetch by journal issn
    * **min_date=(str)**: Fetch paper that appears after min_date (YYYY-MM-DD)
    * **min_date=(str)**: Fetch paper that appears before max_date (YYYY-MM-DD)
    * **time_span=(str)**: Fetch paper that appears in the last time_span days
    * **search=(str)**: Fetch papers on title, journal, authors first and last names

    ** Ordering: **

    * **ordering=[+|-](str) **: order results by field. Choices are: 'date_co' or 'altmetric__score'

    ** Sub-routes: **

    * **[GET] /papers/<id\>/neighbors**: Paper neighbors
    * **[GET] /papers/<id\>/related-threads**: Related threads
        * **time_span=(int)**: Fetch neighbors papers/threads published in the past time_span days (default=60)

    """

    queryset = Paper.objects.filter(is_trusted=True)
    serializer_class = {
        'default': PaperSerializer,
        'nested': PaperNestedSerializer,
    }
    permission_classes = (permissions.AllowAny,
                          )
    filter_backends = (
        DisabledHTMLContextualFilterBackend,
        DisabledHTMLSearchFilterBackend,
        filters.OrderingFilter,
    )
    filter_class = PaperFilter
    ordering_fields = ('date_co',
                       'altmetric__score')
    search_fields = ('title',
                     'journal__title',
                     'authors__first_name',
                     'authors__last_name')

    @method_decorator(cache_page(3600 * 24))
    def list(self, request, *args, **kwargs):
        return super(PaperViewSet, self).list(request, *args, **kwargs)

    @detail_route(methods=['get'])
    def neighbors(self, request, pk=None):
        time_span = int(
            self.request.query_params.get(
                'time_span',
                settings.LIBRARY_DEFAULT_NEIGHBORS_TIMESPAN
            )
        )
        instance = get_object_or_404(self.queryset, pk=pk)
        self.check_object_permissions(request, instance)
        neighbors = instance.get_neighbors(time_span)
        neighbors = self.get_serializer_class()\
            .setup_eager_loading(neighbors)
        serializer = self.get_serializer(neighbors, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @detail_route(methods=['get'], url_path='related-threads')
    def related_threads(self, request, pk=None):
        time_span = int(
            self.request.query_params.get(
                'time_span',
                settings.LIBRARY_DEFAULT_NEIGHBORS_TIMESPAN
            )
        )
        instance = get_object_or_404(self.queryset, pk=pk)
        self.check_object_permissions(request, instance)
        threads = instance.get_related_threads(time_span)
        serializer = ThreadSerializer(
            threads,
            context=self.get_serializer_context(),
            many=True
        )
        return Response(serializer.data, status=status.HTTP_200_OK)

    @detail_route(methods=['get'], url_path='json-ld')
    def json_ld(self, request, pk=None):
        """Provide JSON-LD microdata"""
        instance = get_object_or_404(self.queryset, pk=pk)
        self.check_object_permissions(request, instance)
        serializer = JsonLdSerializer(instance,
                                      context=self.get_serializer_context())
        return Response(serializer.data, status=status.HTTP_200_OK)

    def get_queryset(self):
        queryset = super(PaperViewSet, self).get_queryset()

        # Disabling fetch of entire library with no arguments
        if self.action == 'list' and not self.request.query_params:
            return Paper.objects.none()

        # Altmetric flag to prefetch altmetric score for ordering
        altmetric_flag = False
        if 'altmetric' in self.request.query_params.get('ordering', ''):
            altmetric_flag = True

        return self.get_serializer_class()\
            .setup_eager_loading(queryset,
                                 user=self.request.user,
                                 request=self.request,
                                 altmetric=altmetric_flag)


class MyPaperViewSet(PaperViewSet):
    """
    My Papers

    ### Routes ###

    * **[GET, POST] /my-papers/**: List of papers
    * **[GET] /my-papers/filters**: Filter list for request papers list
    * **[GET, PUT, PATCH] /my-papers/<id\>/**: Paper instance
    * **[GET] /my-papers/<id\>/neighbors**: Paper neighbors
    * **[GET] /my-papers/<id\>/related-threads**: Related threads
    * **[GET] /my-papers/<id\>/json-ld**: JSON-LD based on schema.org ScholarlyArticle item

    ### Additional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested', 'default'

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested', 'default'
    * **doi=(str)**: Fetch by DOI
    * **pmi=(str)**: Fetch by Pub Med Id
    * **arxiv=(str)**: Fetch by Arxiv id
    * **pii=(str)**: Fetch by Elsevier id
    * **title=(str)**: Fetch by title
    * **journal=(str)**: Fetch by journal title
    * **issn=(str)**: Fetch by journal issn
    * **min_date=(str)**: Fetch paper that appears after min_date (YYYY-MM-DD)
    * **min_date=(str)**: Fetch paper that appears before max_date (YYYY-MM-DD)
    * **pinned=(int)**: Fetch only **pinned** (if 1) or **non pinned** (if 0) papers for current user (default = Null)
    * **added=(int)**: Fetch only **added** (if 1) or **non added** (if 0) papers for current user (default = Null)
    * **trashed=(int)**: Fetch only **trashed** (if 1) or **non trashed** (if 0) papers for current user (default = Null)
    * **banned=(int)**: Fetch only **banned** (if 1) or **non banned** (if 0) papers for current user (default = Null)
    * **scored=(int)&type=(str)&feed=(str)**: Fetch only **scored** (if 1) or **non scored** (if 0) papers for current
    user, type ['stream'|'trend'] and feed name (default = (Null, 'stream', 'main'))
    * **time_span=(int)**: Fetch only papers published in the past time_span days
    * **journal_id=(int)**: Filter papers by journal_id
    * **author_id=(int)**: Filter papers by author_id
    * **search=(str)**: Filter papers on title, journal, authors first and last names
    """

    queryset = Paper.objects.all()
    serializer_class = {
        'default': PaperSerializer,
        'nested': PaperNestedSerializer
    }
    exclude_action_serializers = {
        # 'list': ['nested'],
    }
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
                          )

    filter_backends = (DisabledHTMLContextualFilterBackend,
                       DisabledHTMLSearchFilterBackend,
                       )
    filter_class = MyPaperFilter
    search_fields = ('title',
                     'journal__title',
                     'authors__first_name',
                     'authors__last_name')

    SIZE_MAX_JOURNAL_FILTER = 40
    SIZE_MAX_AUTHOR_FILTER = 40
    NEIGHBORS_TIME_SPAN = 60
    TIME_SPAN_DEFAULT = 30
    AUTHOR_COUNT_BUMPER = 1.  # in percent of all journals in user fingerprint
    JOURNAL_COUNT_BUMPER = 5.  # in percent of all journals in user fingerprint

    @method_decorator(never_cache)
    def list(self, request, *args, **kwargs):
        return super(PaperViewSet, self).list(request, *args, **kwargs)

    def get_paper_id(self):
        return self.kwargs['pk']

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def store_controls(self, request):
        # Store persistent user control states
        search = request.query_params.get('search', None)
        scored = request.query_params.get('scored', None)
        time_span = int(request.query_params.get(
            'time_span',
            self.TIME_SPAN_DEFAULT)
        )
        pin = request.query_params.get('pinned', None)
        if scored == '1':
            request.session['feeds-control-states'] = {
                'time_span': time_span,
                # 'search': search,
                'pin': 1 if pin == '1' else 0
            }
        else:
            request.session['library-control-states'] = {
                'time_span': time_span,
                # 'search': search,
                'pin': 1 if pin == '1' else 0
            }

    def get_queryset(self):
        # Store session control
        self.store_controls(self.request)

        queryset = super(PaperViewSet, self).get_queryset()
        return self.get_serializer_class()\
            .setup_eager_loading(queryset,
                                 user=self.request.user,
                                 request=self.request)

    # @method_decorator(cache_page(60 * 60))
    @list_route(methods=['get'])
    def filters(self, request):
        queryset = super(PaperViewSet, self).get_queryset()
        queryset = queryset.select_related('journal')
        queryset = queryset.prefetch_related('authors')

        if self.request.query_params.get('scored') == 1:
            data = list(self.filter_queryset(queryset)[:settings.FEED_N_FIRST_PAPERS_ONLY])
        else:
            data = list(self.filter_queryset(queryset))

        authors = []
        journals = []
        if data and request.user.is_authenticated() and request.user.type == USER_INDIVIDUAL:

            # Journals
            js = [d.journal for d in data if d.journal and d.journal.title]
            js_count = Counter(js).most_common()
            for j, c in js_count[:self.SIZE_MAX_JOURNAL_FILTER]:
                j.count = c
                journals.append(j)

            # Authors
            da = [auth for d in data for auth in d.authors.all()
                  if not auth.first_name.strip() == '' or
                  not auth.last_name.strip() == '']
            as_count = Counter(da).most_common()

            for a, c in as_count[:self.SIZE_MAX_AUTHOR_FILTER]:
                a.count = c
                authors.append(a)

            if self.request.query_params.get('scored') == 1:

                # Get user fingerprint
                f = request.user.fingerprint\
                    .values('authors_ids', 'authors_counts', 'journals_ids', 'journals_counts')[0]
                # Get journal_id and author_id that are bumping up the list
                sac = sum(f['authors_counts']) or 1
                sjc = sum(f['journals_counts']) or 1
                aids_bumping = [aid for i, aid in enumerate(f['authors_ids'])
                                if f['authors_counts'][i]/sac > self.AUTHOR_COUNT_BUMPER/100.]
                jids_bumping = [jid for i, jid in enumerate(f['journals_ids'])
                                if f['journals_counts'][i]/sjc > self.JOURNAL_COUNT_BUMPER/100.]

                # Bump journals
                for i, j in enumerate(journals):
                    if j.id in jids_bumping:
                        journals.insert(0, journals.pop(i))

                # Bump authors
                for i, a in enumerate(authors):
                    if a.id in aids_bumping:
                        authors.insert(0, authors.pop(i))
        else:
            authors = []
            journals = []

        kwargs = {'context': self.get_serializer_context()}
        serializer = PaperFilterSerializer({'authors': authors,
                                            'journals': journals}, **kwargs)

        return Response(serializer.data,  status=status.HTTP_200_OK)


class JournalViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Paper

    ### Routes ###

    * **[GET] /journals/**: List of journals
    * **[GET] /journals/<id\>/**: Journal instance
    """
    queryset = Journal.objects.all()
    serializer_class = JournalSerializer
    permission_classes = (IsReadOnlyRequest, )


class AuthorViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Author

    ### Routes ###

    * [GET] /authors/: List of authors
    * [GET] /authors/<id\>/: Author instance
    """
    queryset = Author.objects.all()
    serializer_class = AuthorSerializer
    permission_classes = (IsReadOnlyRequest, )


class PaperStateViewSet(MultiSerializerMixin,
                        mixins.CreateModelMixin,
                        mixins.ListModelMixin,
                        mixins.UpdateModelMixin,
                        mixins.RetrieveModelMixin,
                        viewsets.GenericViewSet):
    """
    Paper User state: Relation between a [User][ref1] and a [Paper][ref2]

    ### Routes ###

    * **[GET, POST] /states/**: List of Paper/User states
    * **[GET, PUT, PATCH] /states/<id\>/**: Paper/User state instance
    * **[DELETE] /states/empty-trash**: Empty trash

    [ref1]: /api/v1/user/users/
    [ref2]: /api/v1/library/papers/

    """

    queryset = PaperUser.objects.all()
    serializer_class = {
        'default': PaperUserSerializer,
    }
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
                          IsOwner,
                          )

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return PaperUser.objects.filter(user=self.request.user)
        return PaperUser.objects.all()

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    @list_route(methods=['delete'], url_path='empty-trash')
    def empty_trash(self, request):
        from etalia.users.models import UserLibPaper
        pu_qs = PaperUser.objects.filter(
            user=request.user.id,
            store=PAPER_TRASHED
        )
        pids = [d.paper_id for d in pu_qs]
        ulp_qs = UserLibPaper.objects.filter(
            userlib__user=request.user.id,
            paper_id__in=pids
        )
        pu_qs.delete()
        ulp_qs.delete()
        return Response('Trash empty', status=status.HTTP_200_OK)
