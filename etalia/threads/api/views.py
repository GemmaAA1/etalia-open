# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import operator
from functools import reduce
from collections import Counter
from django.contrib.auth import get_user_model

from rest_framework import viewsets, permissions, mixins, status
from rest_framework.exceptions import ParseError
from rest_framework.response import Response
from rest_framework.decorators import detail_route, list_route

from django.db.models import Q
from django.utils import timezone

from etalia.core.api.permissions import IsThreadMember, IsOwner, \
    IsOwnerOrReadOnly, IsNOTThreadMember, ThreadIsNotYetPublished, \
    ThreadIsPublished, ThreadIsNotYetPublishedIsOwnerIfDeleteMethod, \
    IsToUserOrOwnersReadOnly
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser, ThreadUserInvite
from ..constant import THREAD_JOINED, THREAD_LEFT, THREAD_PINNED, THREAD_BANNED, \
    THREAD_PRIVACIES, THREAD_PUBLIC, THREAD_PRIVATE, THREAD_INVITE_PENDING, \
    THREAD_INVITE_DECLINED, THREAD_INVITE_ACCEPTED
from .serializers import \
    ThreadPostSerializer, ThreadCommentSerializer, ThreadSerializer, \
    ThreadUserSerializer, ThreadNestedSerializer, ThreadPostNestedSerializer, \
    ThreadFilterSerializer, ThreadUserInviteSerializer, \
    ThreadUserInviteUpdateSerializer, ThreadUserUpdateSerializer
from etalia.core.api.mixins import MultiSerializerMixin

User = get_user_model()


class ThreadViewSet(MultiSerializerMixin,
                    viewsets.ModelViewSet):

    """
    Threads

    ### Routes ###

    * **[GET, POST] /threads/**: List of threads
    * **[GET] /threads/filters**: Filter list for request threads list
    * **[GET, PUT, PATCH] /threads/<id\>/**: Thread instance
    * **[PUT, PATCH] /threads/<id\>/publish**: Publish Thread
    * **[GET] /threads/<id\>/neighbors**: Thread neighbors

    ### Optional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',
    * **owned=(int)**: Fetch only **owned** (if 1) or **non owned** (if 0) threads for current user (default = Null)
    * **pinned=(int)**: Fetch only **pinned** (if 1) or **non pinned** (if 0) threads for current user (default = Null)
    * **joined=(int)**: Fetch only **joined** (if 1) or **non joined** (if 0) threads for current user (default = Null)
    * **left=(int)**: Fetch only **left** (if 1) or **non left** (if 0) threads for current user (default = Null)
    * **banned=(int)**: Fetch only **banned** (if 1) or **non banned** (if 0) threads for current user (default = Null)
    * **published=(int)**: Fetch only **published** (if 1) or **non published** (if 0) threads for current user (default = Null)
    * **scored=(int)&feed=(str)**: Fetch only **scored** (if 1) or **non scored** (if 0) threads for current user and specific feed (default = (Null, 'main')
    * **private=(int)**: Fetch only **private** (if 1) or **non private** (if 0) threads for current user (default = Null)    
    * **public=(int)**: Fetch only **public** (if 1) or **non public** (if 0) threads (default = Null)
    * **invited=(int)**: Fetch only **invited** (if 1) or **non invited** (if 0) threads for current user (default = Null)
    * **invited-pending=(int)**: Fetch only **pending invites** (if 1) or **non invites** (if 0) threads for current user (default = Null)
    * **invited-accepted=(int)**: Fetch only **accepted invite** (if 1) or **declined** (if 0) threads for current user (default = Null)
    * **time-span=(int)**: Fetch only threads published in the past time-span days
    * **sort-by=(str)**: Sort threads by: 'date', 'score', 'published-date' (default = published-date)
    * **user_id[]=(int)**: Filter threads by user_id
    * **type[]=(int)**: Filter threads by type
    * **search=(str)**: Filter threads on title, owner first and last names

    ** Sub-routes: **

    * **[GET] /threads/<id\>/neighbors**: Thread neighbors:

        * **time-span=(int)**: Fetch neighbors threads published in the past time-span days (default=60)
    """

    queryset = Thread.objects.all()
    serializer_class = {
        'default': ThreadSerializer,
        'nested': ThreadNestedSerializer,
    }
    exclude_action_serializers = {
        'list': ['nested'],
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsOwnerOrReadOnly,
                          ThreadIsNotYetPublishedIsOwnerIfDeleteMethod)

    query_params_props = {
        'owned': {'type': int, 'min': 0, 'max': 1},
        'joined': {'type': int, 'min': 0, 'max': 1},
        'pinned': {'type': int, 'min': 0, 'max': 1},
        'left': {'type': int, 'min': 0, 'max': 1},
        'banned': {'type': int, 'min': 0, 'max': 1},
        'published': {'type': int, 'min': 0, 'max': 1},
        'scored': {'type': int, 'min': 0, 'max': 1},
        'private': {'type': int, 'min': 0, 'max': 1},
        'public': {'type': int, 'min': 0, 'max': 1},
        'invited': {'type': int, 'min': 0, 'max': 1},
        'invited-pending': {'type': int, 'min': 0, 'max': 1},
        'invited-accepted': {'type': int, 'min': 0, 'max': 1},
        'time-span': {'type': int, 'min': -1, 'max': 1e6},
        'view': {'type': str},
        'sort-by': {'type': str},
        'type[]': {'type': list},
        'user_id[]': {'type': list},
    }

    size_max_user_filter = 40
    neighbors_time_span = 60

    def get_thread_id(self):
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
        feed_name = self.request.query_params.get('feed', 'main')
        bool_filters_def = {
            'owned': {
                'toggle': Q(user=self.request.user)
            },
            'joined': {
                'base': None,
                'toggle': Q(threaduser__user=self.request.user) &
                          Q(threaduser__participate=THREAD_JOINED)
            },
            'left': {
                'base': None,
                'toggle': Q(threaduser__user=self.request.user) &
                          Q(threaduser__participate=THREAD_LEFT)
            },
            'pinned': {
                'base': None,
                'toggle': Q(threaduser__user=self.request.user) &
                          Q(threaduser__watch=THREAD_PINNED)
            },
            'banned': {
                'base': None,
                'toggle': Q(threaduser__user=self.request.user) &
                          Q(threaduser__watch=THREAD_BANNED)
            },
            'published': {
                'base': Q(user=self.request.user),
                'toggle': ~Q(published_at=None)
            },
            'scored': {
                'toggle': Q(threadfeedthreads__threadfeed__name=feed_name) &
                          Q(threadfeedthreads__threadfeed__user=self.request.user)
            },
            'private': {
                'toggle': Q(privacy=THREAD_PRIVATE)
            },
            'public': {
                'toggle': Q(privacy=THREAD_PUBLIC)
            },
            'invited': {
                'toggle': Q(threaduserinvite__to_user=self.request.user)
            },
            'invited-pending': {
                'base': None,
                'toggle': Q(threaduserinvite__to_user=self.request.user) &
                          Q(threaduserinvite__status=THREAD_INVITE_PENDING)
            },
            'invited-accepted': {
                'base': None,
                'toggle': Q(threaduserinvite__to_user=self.request.user) &
                          Q(threaduserinvite__status=THREAD_INVITE_ACCEPTED)
            },
        }

        # ordering mapping
        order_by_map = {
            'published-date': '-published_at',
            'date': '-threaduser__modified',
            'score': '-threadscore__score',
        }

        # Baseline query
        query_args = [
            ~(Q(published_at=None) & ~Q(user=self.request.user)),
            ~(Q(privacy=THREAD_PRIVATE) &
              ~(Q(threaduser__user=self.request.user) & Q(threaduser__participate=THREAD_JOINED)))
        ]

        # boolean filters
        for key, props in bool_filters_def.items():
            param = self.request.query_params.get(key, 'null')
            if not param == 'null':
                if props.get('base', None):
                    query_args.append(props['base'])
                    # base = reduce(operator.and_, props['base'])
                    # queryset = queryset.filter(base)
                if props.get('toggle', None):
                    # toggle = reduce(operator.and_, props['toggle'])
                    if param == '1':
                        query_args.append(props['toggle'])
                    #     queryset = queryset.filter(toggle)
                    elif param == '0':
                        query_args.append(operator.not_(props['toggle']))
                        # queryset = queryset.exclude(toggle)

        # Thread Types
        thread_types = [int(id_) for id_ in self.request.query_params.getlist('type[]', None)]
        if thread_types:
            query_args.append(Q(type__in=thread_types))

        # Owner of thread
        uids = [int(id_) for id_ in self.request.query_params.getlist('user_id[]', None)]
        if uids:
            query_args.append(
                Q(user_id__in=uids)
                # | Q(posts__user_id__in=uids) \
                # | Q(posts__comments__user_id__in=uids)
            )

        # time-span filter
        time_span = self.request.query_params.get('time-span', 'null')
        if not time_span == 'null':
            cutoff_datetime = timezone.now() - timezone.timedelta(
                days=int(time_span))
            query_args.append(Q(published_at__gt=cutoff_datetime))

        # search
        search = self.request.query_params.get('search', 'null')
        if not search == 'null':
            subset = []
            for word in search.split():
                subset.append(Q(title__icontains=word) |
                              Q(user__first_name__icontains=word) |
                              Q(user__last_name__icontains=word))
            if subset:
                query_args.append(reduce(operator.and_, subset))

        queryset = Thread.objects.all().filter(reduce(operator.and_, query_args)).distinct()

        order_by = self.request.query_params.get('sort-by', 'null')
        if not order_by == 'null':
            queryset = queryset.order_by(order_by_map.get(order_by))

        return queryset

    @detail_route(methods=['put', 'patch'],
                  permission_classes=(ThreadIsNotYetPublished,))
    def publish(self, request, pk=None):
        instance = self.get_object()
        self.check_object_permissions(request, instance)
        instance.publish()
        serializer = self.get_serializer(instance)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @detail_route(methods=['get'],
                  permission_classes=(ThreadIsPublished, ))
    def neighbors(self, request, pk=None):
        time_span = int(self.request.query_params.get('time-span',
                                                      self.neighbors_time_span))
        instance = Thread.objects.get(id=pk)
        self.check_object_permissions(request, instance)
        neighbors = instance.get_neighbors(time_span)
        serializer = self.get_serializer(neighbors, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @list_route(methods=['get'])
    def filters(self, request):
        tids = self.get_queryset().values_list('id', flat=True)
        values = ', '.join(['({0})'.format(i) for i in tids])
        du = []
        if values:
            qu = User.objects.raw(
                        "SELECT * "
                        "FROM users_user u "
                        "LEFT JOIN threads_thread t ON u.id = t.user_id "
                        "WHERE t.id IN (VALUES {0}) ".format(values))
            du = list(qu)

        us_count = Counter(du).most_common()
        users = []
        for u, c in us_count[:self.size_max_user_filter]:
            u.count = c
            users.append(u)

        kwargs = {'context': self.get_serializer_context()}
        serializer = ThreadFilterSerializer({'users': users}, **kwargs)
        return Response(serializer.data,  status=status.HTTP_200_OK)


class ThreadPostViewSet(MultiSerializerMixin,
                        viewsets.ModelViewSet):
    """
    Thread Post: Post on thread

    ### Routes ###

    * **[GET, POST] /posts/**: List of posts
    * **[GET, PUT, PATCH] /posts/<id\>/**: Post instance

    ### Optional Kwargs ###

    ** Detail:**

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',
    * **thread_id=(int)**: Filter comments related to Thread thread_id

    """
    queryset = ThreadPost.objects.all()
    serializer_class = {
        'default': ThreadPostSerializer,
        'nested': ThreadPostNestedSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsOwnerOrReadOnly,
                          IsThreadMember
                          )

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def get_queryset(self):
        if self.action == 'list':
            threads_joined = ThreadUser.objects \
                .filter(user=self.request.user, participate=True) \
                .values('thread')
            queryset = ThreadPost.objects.filter(thread__in=threads_joined)

            # filter based on thread_id
            thread_id = self.request.query_params.get('thread_id', 'null')
            if not thread_id == 'null':
                queryset = queryset.filter(thread_id=thread_id)

            return queryset
        return ThreadPost.objects.all()


class ThreadCommentViewSet(MultiSerializerMixin,
                           viewsets.ModelViewSet):
    """
    Thread Comment: Comments on thread post

    ### Routes ###

    * **[GET, POST] /comments/**: List of comments
    * **[GET, PUT, PATCH] /comments/<id\>/**: Comment instance

    ### Optional Kwargs ###

    ** List: **

    * **post_id=(int)**: Filter comments related to Post post_id

    """

    queryset = ThreadComment.objects.filter()
    serializer_class = {
        'default': ThreadCommentSerializer,
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsOwnerOrReadOnly,
                          IsThreadMember
                          )

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def get_queryset(self):
        if self.action == 'list':
            threads_joined = ThreadUser.objects \
                .filter(user=self.request.user, participate=THREAD_JOINED) \
                .values('thread')

            queryset = ThreadComment.objects.filter(
                post__thread__in=threads_joined)

            # filter based on post_id
            post_id = self.request.query_params.get('post_id', 'null')
            if not post_id == 'null':
                queryset = queryset.filter(post_id=post_id)

            return queryset
        return ThreadComment.objects.all()


class ThreadUserViewSet(MultiSerializerMixin,
                        mixins.CreateModelMixin,
                        mixins.ListModelMixin,
                        mixins.UpdateModelMixin,
                        mixins.RetrieveModelMixin,
                        viewsets.GenericViewSet):
    """
    Thread User state: Relation between a [User][ref1] and a [Thread][ref2]

    ### Routes ###

    * **[GET, POST] /states/**: List of Thread User states
    * **[GET, PUT, PATCH] /states/<id\>/**: Thread User state instance

    [ref1]: /api/v1/user/users/
    [ref2]: /api/v1/thread/threads/

    """

    queryset = ThreadUser.objects.all()
    serializer_class = {
        'default': ThreadUserSerializer,
        'update': ThreadUserUpdateSerializer,
        'partial_update': ThreadUserUpdateSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsOwner,
                          )

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return ThreadUser.objects.filter(user=self.request.user)
        return ThreadUser.objects.all()

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)


class ThreadUserInviteViewSet(MultiSerializerMixin,
                              mixins.CreateModelMixin,
                              mixins.ListModelMixin,
                              mixins.UpdateModelMixin,
                              mixins.RetrieveModelMixin,
                              viewsets.GenericViewSet):

    """
    Thread from_user-to-to_user Invite

    ### Routes ###

    * **[GET, POST] /invites/**: List of Invite
    * **[GET, PUT, PATCH] /invites/<id\>/**: Invite

    ### Optional Kwargs ###

    ** List: **

    * **from-user=(int)**: Filter ThreadUserInvite with from-user=user_id
    * **to-user=(int)**: Filter ThreadUserInvite with to-user=user_id
    * **status=(int)**: Filter ThreadUserInvite with specific status

    """

    queryset = ThreadUserInvite.objects.all()
    serializer_class = {
        'default': ThreadUserInviteSerializer,
        'update': ThreadUserInviteUpdateSerializer,
        'partial_update': ThreadUserInviteUpdateSerializer,
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsToUserOrOwnersReadOnly,
                          )

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            queryset = ThreadUserInvite.objects.filter(
                Q(from_user=self.request.user) |
                Q(to_user=self.request.user))
            # filter from_user
            from_user = self.request.query_params.get('from-user', 'null')
            if not from_user == 'null':
                queryset = queryset.filter(from_user=from_user)
            # filter to_user
            to_user = self.request.query_params.get('to-user', 'null')
            if not to_user == 'null':
                queryset = queryset.filter(to_user=to_user)
            # filter status
            status = self.request.query_params.get('status', 'null')
            if not status == 'null':
                queryset = queryset.filter(status=status)

            return queryset

        return ThreadUserInvite.objects.all()

    def perform_create(self, serializer):
        serializer.save(from_user=self.request.user)
