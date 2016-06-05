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
from django.db import connection

from etalia.core.api.permissions import IsReadOnlyRequest
from etalia.core.api.mixins import MultiSerializerMixin
from etalia.core.api.permissions import IsOwner

from ..models import Paper, Author, Journal, PaperUser
from etalia.users.models import UserLibPaper
from ..constants import PAPER_BANNED, PAPER_PINNED

from .serializers import PaperSerializer, JournalSerializer, AuthorSerializer, \
    PaperNestedSerializer, PaperFilterSerializer, PaperStateSerializer


class PaperViewSet(MultiSerializerMixin, viewsets.ModelViewSet):
    """
    Paper

    Papers

    ### Routes ###

    * **[GET, POST] /papers/**: List of papers
    * **[GET] /papers/filters**: Filter list for request papers list
    * **[GET, PUT, PATCH] /papers/<id\>/**: Paper instance
    * **[GET] /papers/<id\>/neighbors**: Paper neighbors

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

    ** Sub-routes: **

    * **[GET] /papers/<id\>/neighbors**: Paper neighbors:

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

    size_max_journal_filter = 40
    size_max_author_filter = 40
    neighbors_time_span = 60

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
                'query': [Q(userlib_paper__user=self.request.user)],
            },
            'trashed': {
                'query': [Q(userlib_paper__user=self.request.user),
                          Q(userlib_paper__is_trashed=True)],
            },
            'pinned': {
                'query': [Q(paperuser__user=self.request.user),
                          Q(paperuser__watch=PAPER_PINNED)],
            },
            'banned': {
                'query': [Q(paperuser__user=self.request.user),
                          Q(paperuser__watch=PAPER_BANNED)],
            },
        }

        # base queryset
        queryset = self.queryset

        # boolean filters
        for key, props in bool_filters_def.items():
            param = self.request.query_params.get(key, None)
            query = reduce(operator.and_, props['query'])
            if param == '1':
                queryset = queryset.filter(query)
            elif param == '0':
                queryset = queryset.exclude(query)

        # Filters
        jids = [int(id_) for id_ in self.request.query_params.getlist('journal_id[]', None)]
        if jids:
            queryset = queryset.filter(
                Q(journal_id__in=jids)
            )
        aids = [int(id_) for id_ in self.request.query_params.getlist('author_id[]', None)]
        if aids:
            queryset = queryset.filter(
                Q(authors__in=aids)
            )

        # Paper feeds
        scored = self.request.query_params.get('scored', None)
        type = self.request.query_params.get('type', 'stream')
        feed_name = self.request.query_params.get('feed', 'main')
        time_span = self.request.query_params.get('time-span', None)
        if scored == '1':
            if type == 'stream':
                queryset = queryset.filter(
                    Q(streampapers__stream__name=feed_name),
                    Q(streampapers__stream__user=self.request.user))\
                    .order_by('-streampapers__score')
            elif type == 'trend':
                queryset = queryset.filter(
                    Q(trendpapers__trend__name=feed_name),
                    Q(trendpapers__trend__user=self.request.user))\
                    .order_by('-trendpapers__score')

            # time-span filter
            if time_span:
                cutoff_datetime = timezone.now() - timezone.timedelta(
                    days=int(time_span))
                if type == 'stream':
                    queryset = queryset.filter(streampapers__date__gt=cutoff_datetime)
                elif type == 'trend':
                    queryset = queryset.filter(trendpapers__date__gt=cutoff_datetime)

            return queryset

        queryset = queryset.order_by('-userlib_paper__date_created')

        return queryset

    @detail_route(methods=['get'])
    def neighbors(self, request, pk=None):
        time_span = int(self.request.query_params.get('time-span',
                                                      self.neighbors_time_span))
        instance = self.get_object()
        self.check_object_permissions(request, instance)
        neighbors = instance.get_neighbors(time_span)
        serializer = self.get_serializer(neighbors, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @list_route(methods=['get'])
    def filters(self, request):
        queryset = self.get_queryset()

        data = queryset.select_related('journal')
        pids = [d.id for d in data]

        # Journals
        js = [d.journal for d in data]
        js_count = Counter(js).most_common()
        journals = []
        for j, c in js_count[:self.size_max_journal_filter]:
            j.count = c
            journals.append(j)

        # Authors
        values = ', '.join(['({0})'.format(i) for i in pids])
        qa = Author.objects.raw(
                    "SELECT * "
                    "FROM library_author a "
                    "LEFT JOIN library_authorpaper ap ON a.id = ap.author_id "
                    "WHERE ap.paper_id IN (VALUES {0}) ".format(values))
        da = list(qa)
        as_count = Counter(da).most_common()
        authors = []
        for a, c in as_count[:self.size_max_author_filter]:
            a.count = c
            authors.append(a)

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

    [ref1]: /api/v1/user/users/
    [ref2]: /api/v1/library/papers/

    """

    queryset = PaperUser.objects.all()
    serializer_class = {
        'default': PaperStateSerializer,
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