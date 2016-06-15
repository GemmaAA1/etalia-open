# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import operator
from collections import Counter
from functools import reduce

from rest_framework import viewsets, permissions, status, mixins
from rest_framework.exceptions import ParseError
from rest_framework.decorators import detail_route, list_route
from rest_framework.response import Response

from django.db.models import Q
from django.utils import timezone

from etalia.core.api.permissions import IsReadOnlyRequest
from etalia.core.api.mixins import MultiSerializerMixin
from etalia.core.api.permissions import IsOwner
from etalia.threads.api.serializers import ThreadSerializer

from ..models import Paper, Author, Journal, PaperUser
from ..constants import PAPER_BANNED, PAPER_PINNED, PAPER_ADDED, PAPER_TRASHED

from .serializers import PaperSerializer, JournalSerializer, AuthorSerializer, \
    PaperNestedSerializer, PaperFilterSerializer, PaperUserSerializer, \
    PaperUserUpdateSerializer


class PaperViewSet(MultiSerializerMixin,
                   viewsets.ReadOnlyModelViewSet):
    """
    Paper

    Papers

    ### Routes ###

    * **[GET, POST] /papers/**: List of papers
    * **[GET] /papers/filters**: Filter list for request papers list
    * **[GET, PUT, PATCH] /papers/<id\>/**: Paper instance
    * **[GET] /papers/<id\>/neighbors**: Paper neighbors
    * **[GET] /papers/<id\>/related-threads**: Related threads

    ### Optional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',
    * **pinned=(int)**: Fetch only **pinned** (if 1) or **non pinned** (if 0) papers for current user (default = Null)
    * **added=(int)**: Fetch only **added** (if 1) or **non added** (if 0) papers for current user (default = Null)
    * **trashed=(int)**: Fetch only **trashed** (if 1) or **non trashed** (if 0) papers for current user (default = Null)
    * **banned=(int)**: Fetch only **banned** (if 1) or **non banned** (if 0) papers for current user (default = Null)
    * **scored=(int)&type=(str)&feed=(str)**: Fetch only **scored** (if 1) or **non scored** (if 0) papers for current
    user, type ['stream'|'trend'] and feed name (default = (Null, 'stream', 'main'))
    * **time-span=(int)**: Fetch only papers published in the past time-span days
    * **journal_id[]=(int)**: Filter papers by journal_id
    * **author_id[]=(int)**: Filter papers by author_id
    * **search=(str)**: Filter papers on title, journal, authors first and last names

    ** Sub-routes: **

    * **[GET] /papers/<id\>/neighbors**: Paper neighbors
    * **[GET] /papers/<id\>/related-threads**: Related threads

        * **time-span=(int)**: Fetch neighbors papers published in the past time-span days (default=60)
    """

    queryset = Paper.objects.all()
    serializer_class = {
        'default': PaperSerializer,
        'nested': PaperNestedSerializer
    }
    exclude_action_serializers = {
        # 'list': ['nested'],
    }
    permission_classes = (permissions.IsAuthenticated,
                          )

    query_params_props = {
        'added': {'type': int, 'min': 0, 'max': 1},
        'pinned': {'type': int, 'min': 0, 'max': 1},
        'trashed': {'type': int, 'min': 0, 'max': 1},
        'banned': {'type': int, 'min': 0, 'max': 1},
        'scored': {'type': int, 'min': 0, 'max': 1},
        'time-span': {'type': int, 'min': 0, 'max': 1e6},
        'view': {'type': str},
        'sort-by': {'type': str},
        'journal_id[]': {'type': list},
        'author_id[]': {'type': list},
    }

    SIZE_MAX_JOURNAL_FILTER = 40
    SIZE_MAX_AUTHOR_FILTER = 40
    NEIGHBORS_TIME_SPAN = 60
    TIME_SPAN_DEFAULT = 30
    AUTHOR_COUNT_BUMPER = 1.  # in percent of all journals in user fingerprint
    JOURNAL_COUNT_BUMPER = 5.  # in percent of all journals in user fingerprint

    def get_paper_id(self):
        return self.kwargs['pk']

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def validate_query_params(self):
        for key, props in self.query_params_props.items():
            if props['type'] == list:
                param = self.request.query_params.getlist(key, None)
            else:
                param = self.request.query_params.get(key, None)
            if param:
                try:
                    val = props['type'](param)
                    if 'min' in props:
                        assert val >= props['min'], \
                            'is inferior to allowed min ({min})'.format(
                                min=props['min'])
                    if 'max' in props:
                        assert val <= props['max'], \
                            'is superior to allowed max ({max})'.format(
                                max=props['max'])
                except (ValueError, AssertionError) as err:
                    raise ParseError('{key}: {err}'.format(key=key, err=err))

    def get_queryset(self):

        # validate query_params
        self.validate_query_params()

        # Bool filters definition
        bool_filters_def = {
            'added': {
                'base': None,
                'toggle': Q(paperuser__user=self.request.user) &
                          Q(paperuser__store=PAPER_ADDED)
            },
            'trashed': {
                'base': None,
                'toggle': Q(paperuser__user=self.request.user) &
                          Q(paperuser__store=PAPER_TRASHED)
            },
            'pinned': {
                'base': None,
                'toggle': Q(paperuser__user=self.request.user) &
                          Q(paperuser__watch=PAPER_PINNED)
            },
            'banned': {
                'base': None,
                'toggle': Q(paperuser__user=self.request.user) &
                          Q(paperuser__watch=PAPER_BANNED)
            },
        }

        # base queryset
        queryset = self.queryset
        query_args = []
        order_by = '-userlib_paper__date_created'

        # boolean filters
        for key, props in bool_filters_def.items():
            param = self.request.query_params.get(key, 'null')
            if not param == 'null':
                if props.get('base', None):
                    query_args.append(props['base'])
                if props.get('toggle', None):
                    if param == '1':
                        query_args.append(props['toggle'])
                    elif param == '0':
                        query_args.append(~props['toggle'])

        # Filters
        jids = [int(id_) for id_ in self.request.query_params.getlist('journal_id[]', None)]
        if jids:
            query_args.append(Q(journal_id__in=jids))
        aids = [int(id_) for id_ in self.request.query_params.getlist('author_id[]', None)]
        if aids:
            query_args.append(Q(authors__in=aids))

        # search
        search = self.request.query_params.get('search', 'null')
        if not search == 'null':
            subset = []
            for word in search.split():
                subset.append(
                    Q(title__icontains=word) |
                    Q(journal__title__icontains=word) |
                    Q(authors__first_name__icontains=word) |
                    Q(authors__last_name__icontains=word)
                )
            if subset:
                query_args.append(reduce(operator.and_, subset))

        # Paper feeds
        scored = self.request.query_params.get('scored', 'null')
        type = self.request.query_params.get('type', 'stream')
        feed_name = self.request.query_params.get('feed', 'main')
        time_span_str = self.request.query_params.get('time-span',
                                                  self.TIME_SPAN_DEFAULT)
        try:
            time_span = int(time_span_str)
        except ValueError:
            time_span = self.TIME_SPAN_DEFAULT
        if scored == '1':
            if type == 'stream':
                query_args.append(
                    Q(streampapers__stream__name=feed_name) &
                    Q(streampapers__stream__user=self.request.user))
                order_by = '-streampapers__score'
            elif type == 'trend':
                query_args.append(
                    Q(trendpapers__trend__name=feed_name) &
                    Q(trendpapers__trend__user=self.request.user))
                order_by = '-trendpapers__score'

            # time-span filter
            if time_span:
                cutoff_datetime = timezone.now() - timezone.timedelta(
                    days=int(time_span))
                if type == 'stream':
                    query_args.append(Q(streampapers__date__gt=cutoff_datetime))
                elif type == 'trend':
                    query_args.append(Q(trendpapers__date__gt=cutoff_datetime))

        if query_args:
            queryset = queryset.filter(reduce(operator.and_, query_args))\
                .order_by(order_by)
        else:
            queryset = queryset.order_by(order_by)

        # Store persistent user control states
        pin = self.request.query_params.get('pinned', 'null')
        if scored == '1':
            self.request.session['feeds-control-states'] = {
                'time-span': time_span,
                'search': search,
                'pin': 1 if pin == '1' else 0
            }
        else:
            self.request.session['library-control-states'] = {
                'time-span': time_span,
                'search': search,
                'pin': 1 if pin == '1' else 0
            }

        return queryset

    @detail_route(methods=['get'])
    def neighbors(self, request, pk=None):
        time_span = int(self.request.query_params.get('time-span',
                                                      self.NEIGHBORS_TIME_SPAN))
        instance = self.get_object()
        self.check_object_permissions(request, instance)
        neighbors = instance.get_neighbors(time_span)
        serializer = self.get_serializer(neighbors, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @detail_route(methods=['get'], url_path='related-threads')
    def related_threads(self, request, pk=None):
        time_span = int(self.request.query_params.get('time-span',
                                                      self.NEIGHBORS_TIME_SPAN))
        instance = self.get_object()
        self.check_object_permissions(request, instance)
        threads = instance.get_related_threads(request.user.id, time_span)
        serializer = ThreadSerializer(threads,
                                      context=self.get_serializer_context(),
                                      many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @list_route(methods=['get'])
    def filters(self, request):
        queryset = self.get_queryset()

        data = queryset.select_related('journal')
        pids = [d.id for d in data]

        authors = []
        journals = []
        if pids:

            # Get user fingerprint
            f = request.user.fingerprint\
                .values('authors_ids', 'authors_counts', 'journals_ids', 'journals_counts')[0]
            # Get journal_id and author_id that are bumping up the list
            sac = sum(f['authors_counts'])
            sjc = sum(f['journals_counts'])
            aids_bumping = [aid for i, aid in enumerate(f['authors_ids'])
                            if f['authors_counts'][i]/sac > self.AUTHOR_COUNT_BUMPER/100.]
            jids_bumping = [jid for i, jid in enumerate(f['journals_ids'])
                            if f['journals_counts'][i]/sjc > self.AUTHOR_COUNT_BUMPER/100.]

            # Journals
            js = [d.journal for d in data if d.journal.title]
            js_count = Counter(js).most_common()
            for j, c in js_count[:self.SIZE_MAX_JOURNAL_FILTER]:
                j.count = c
                journals.append(j)

            # Bump journals
            for i, j in enumerate(journals):
                if j.id in jids_bumping:
                    journals.insert(0, journals.pop(i))

            # Authors
            values = ', '.join(['({0})'.format(i) for i in pids])
            qa = Author.objects.raw(
                        "SELECT * "
                        "FROM library_author a "
                        "LEFT JOIN library_authorpaper ap ON a.id = ap.author_id "
                        "WHERE ap.paper_id IN (VALUES {0}) "
                        "   AND (a.first_name <> '' OR a.last_name <> '')".format(values))
            da = list(qa)
            as_count = Counter(da).most_common()

            for a, c in as_count[:self.SIZE_MAX_AUTHOR_FILTER]:
                a.count = c
                authors.append(a)

            # Bump authors
            for i, a in enumerate(authors):
                if a.id in aids_bumping:
                    authors.insert(0, authors.pop(i))

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
        # 'update': PaperUserUpdateSerializer,
        # 'partial_update': PaperUserUpdateSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
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
        queryset = PaperUser.objects.raw(
            "SELECT pu.id,"
            "       ulp.id AS ulp_id "
            "FROM library_paperuser pu "
            "LEFT JOIN users_userlibpaper ulp ON pu.paper_id = ulp.paper_id "
            "WHERE pu.store = %s "
            "   AND pu.user_id = %s", (PAPER_TRASHED, request.user.id)
        )
        pu_ids = [d.id for d in queryset]
        ulp_ids = [d.ulp_id for d in queryset]
        PaperUser.objects.filter(id__in=pu_ids).delete()
        UserLibPaper.objects.filter(id__in=ulp_ids).delete()
        return Response('Trash empty', status=status.HTTP_200_OK)

