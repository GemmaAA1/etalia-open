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
    IsOwnerOrReadOnly, IsNOTThreadMember, ThreadIsNotYetPublished
from etalia.users.api.serializers import UserFilterSerializer
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser
from ..constant import THREAD_JOINED, THREAD_LEFT, THREAD_PINNED, THREAD_BANNED, \
    THREAD_PRIVACIES, THREAD_PUBLIC, THREAD_PRIVATE, THREAD_INVITE_PENDING, \
    THREAD_INVITE_DECLINED, THREAD_INVITE_ACCEPTED
from .serializers import \
    ThreadPostSerializer, ThreadCommentSerializer, ThreadSerializer, \
    ThreadUserSerializer, ThreadNestedSerializer, ThreadPostNestedSerializer, \
    ThreadCommentNestedSerializer, ThreadFilterSerializer
from etalia.core.api.mixins import MultiSerializerMixin

User = get_user_model()


class ThreadViewSet(MultiSerializerMixin,
                    mixins.CreateModelMixin,
                    mixins.ListModelMixin,
                    mixins.RetrieveModelMixin,
                    mixins.UpdateModelMixin,
                    viewsets.GenericViewSet):
    """
    Threads

    ### Routes ###

    * **[GET, POST] /threads/**: List of threads
    * **[GET, PUT, PATCH] /threads/<id\>/**: Thread instance
    * **[PUT, PATCH] /threads/<id\>/publish**: Publish Thread
    * **[GET] /threads/filters**: Filter list for request threads list

    ### Optional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',
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
                          IsOwnerOrReadOnly)

    query_params_props = {
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
        'time-span': {'type': int, 'min': 0, 'max': 1e6},
        'view': {'type': str},
        'sort-by': {'type': str}
    }

    size_max_user_filter = 40

    def get_thread_id(self):
        return self.kwargs['pk']

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def validate_query_params(self):
        for key, props in self.query_params_props.items():
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
            'joined': {
                'query': [Q(threaduser__user=self.request.user),
                          Q(threaduser__participate=THREAD_JOINED)],
            },
            'pinned': {
                'query': [Q(threaduser__user=self.request.user),
                          Q(threaduser__watch=THREAD_PINNED)],
            },
            'left': {
                'query': [Q(threaduser__user=self.request.user),
                          Q(threaduser__participate=THREAD_LEFT)],
            },
            'banned': {
                'query': [Q(threaduser__user=self.request.user),
                          Q(threaduser__watch=THREAD_BANNED)],
            },
            'published': {
                'query': [Q(user=self.request.user),
                          ~Q(published_at=None)],
            },
            'scored': {
                'query': [Q(threadscore__thread_feed__name=feed_name),
                          Q(threadscore__thread_feed__user=self.request.user)],
            },
            'private': {
                'query': [Q(privacy=THREAD_PRIVATE)],
            },
            'public': {
                'query': [Q(privacy=THREAD_PUBLIC)],
            },
            'invited': {
                'query': [Q(threaduserinvite__to_user=self.request.user)],
            },
            'invited-pending': {
                'query': [Q(threaduserinvite__to_user=self.request.user),
                          Q(threaduserinvite__status=THREAD_INVITE_PENDING)],
            },
            'invited-accepted': {
                'query': [Q(threaduserinvite__to_user=self.request.user),
                          Q(threaduserinvite__status=THREAD_INVITE_ACCEPTED)],
            },
        }

        # ordering mapping
        order_by_map = {
            'published-date': '-published_at',
            'date': '-threaduser__modified',
            'score': '-threadscore__score',
        }

        # base queryset
        queryset = Thread.objects.all() \
            .exclude(Q(published_at=None) & ~Q(user=self.request.user)) \
            .exclude(Q(privacy=THREAD_PRIVATE) & ~(
            Q(threaduser__user=self.request.user) & Q(
                threaduser__participate=THREAD_JOINED)))

        # boolean filters
        for key, props in bool_filters_def.items():
            param = self.request.query_params.get(key, None)
            query = reduce(operator.and_, props['query'])
            if param == '1':
                queryset = queryset.filter(query)
            elif param == '0':
                queryset = queryset.exclude(query)

        # time-span filter
        time_span = self.request.query_params.get('time-span', None)
        if time_span:
            cutoff_datetime = timezone.now() - timezone.timedelta(
                days=int(time_span))
            queryset = queryset.filter(published_at__gt=cutoff_datetime)

        order_by = self.request.query_params.get('sort-by', None)
        if order_by:
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

    @list_route(methods=['get'])
    def filters(self, request):
        queryset = self.get_queryset().select_related('posts__user',
                                                      'posts__comments__user')
        ids = []
        ids += queryset.values_list('posts__user', flat=True)
        ids += queryset.values_list('posts__comments__user', flat=True)
        ids = [id for id in ids if id is not None]
        ids_count = Counter(ids)
        ids_count_top = dict([(pk, count) for pk, count in ids_count.items()][:self.size_max_user_filter])

        clauses = ' '.join(['WHEN id=%s THEN %s' %
                            (pk, count) for pk, count in enumerate(ids_count_top)])
        ordering = 'CASE %s END' % clauses
        users_ordered = list(User.objects\
            .filter(pk__in=ids_count.keys())\
            .extra(select={'ordering': ordering},order_by=('ordering',)))
        for user in users_ordered:
            user.count = ids_count[user.id]
        kwargs = {'context': self.get_serializer_context()}
        serializer = ThreadFilterSerializer({'users': users_ordered}, **kwargs)
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
            thread_id = self.request.query_params.get('thread_id', None)
            if thread_id is not None:
                queryset = queryset.filter(thread_id=thread_id)

            return queryset
        return ThreadPost.objects.all()


class ThreadCommentViewSet(MultiSerializerMixin, viewsets.ModelViewSet):
    """
    Thread Comment: Comments on thread post

    ### Routes ###

    * **[GET, POST] /comments/**: List of comments
    * **[GET, PUT, PATCH] /comments/<id\>/**: Comment instance

    ### Optional Kwargs ###

    ** Detail:**

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',
    * **post_id=(int)**: Filter comments related to Post post_id

    """

    queryset = ThreadComment.objects.filter()
    serializer_class = {
        'default': ThreadCommentSerializer,
        'nested': ThreadCommentNestedSerializer
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
            post_id = self.request.query_params.get('post_id', None)
            if post_id is not None:
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

    **Deprecated**:

    * **[PUT, PATCH] /states/<id\>/join**: User **join** thread
    * **[PUT, PATCH] /states/<id\>/leave**: User **leave** thread
    * **[PUT, PATCH] /states/<id\>/pin**: User **pin** thread
    * **[PUT, PATCH] /states/<id\>/ban**: User **ban** thread

    [ref1]: /api/v1/user/users/
    [ref2]: /api/v1/thread/threads/

    """

    queryset = ThreadUser.objects.filter()
    serializer_class = {
        'default': ThreadUserSerializer,
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

    @detail_route(methods=['put', 'patch'],
                  permission_classes=(IsNOTThreadMember,))
    def join(self, request, pk=None):
        return self.perform_threaduser_action(request, 'join')

    @detail_route(methods=['put', 'patch'],
                  permission_classes=(IsThreadMember,))
    def leave(self, request, pk=None):
        return self.perform_threaduser_action(request, 'leave')

    @detail_route(methods=['put', 'patch'])
    def pin(self, request, pk=None):
        return self.perform_threaduser_action(request, 'pin')

    @detail_route(methods=['put', 'patch'])
    def ban(self, request, pk=None):
        return self.perform_threaduser_action(request, 'ban')

    def perform_threaduser_action(self, request, action):
        """Perform ThreadUser related action (join, pin, etc.)"""
        instance = self.get_object()
        self.check_object_permissions(request, instance)
        method = getattr(instance, action)
        method()
        serializer = self.get_serializer(instance)
        return Response(serializer.data, status=status.HTTP_200_OK)

        # def partial_update(self, request, *args, **kwargs):
        #     kwargs['context'] = self.get_serializer_context()
        #     try:  # pop pk, unrelevant for patch
        #         kwargs.pop('pk')
        #     except KeyError:
        #         pass
        #     serializer_class = self.serializer_class['patch']
        #     serializer = serializer_class(data=request.data, **kwargs)
        #     serializer.is_valid(raise_exception=True)
        #     instance = self.get_object()
        #     self.perform_patch_update(serializer, instance)
        #     serializer = self.get_serializer(instance)
        #     return Response(serializer.data)
