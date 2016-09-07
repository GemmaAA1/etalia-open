# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from collections import Counter

from django.contrib.auth import get_user_model
from django.views.decorators.cache import never_cache
from django.utils.decorators import method_decorator
from django.db.models import Q
from django.shortcuts import get_object_or_404
from django.conf import settings
from rest_framework import viewsets, permissions, mixins, status, filters
from rest_framework.response import Response
from rest_framework.decorators import detail_route, list_route

from etalia.core.api.filters import DisabledHTMLContextualFilterBackend, \
    DisabledHTMLSearchFilterBackend
from etalia.core.api.permissions import IsThreadMember, IsOwner, \
    IsOwnerOrReadOnly, ThreadIsNotYetPublished, \
    ThreadIsPublished, ThreadIsNotYetPublishedIsOwnerIfDeleteMethod, \
    IsToUserOrOwnersReadOnly, IsSessionAuthenticatedOrReadOnly
from etalia.core.api.mixins import MultiSerializerMixin
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser, ThreadUserInvite
from ..constant import THREAD_JOINED, THREAD_LEFT, THREAD_PINNED, THREAD_BANNED, \
    THREAD_PUBLIC, THREAD_PRIVATE, THREAD_INVITE_PENDING, \
    THREAD_INVITE_ACCEPTED
from .serializers import \
    ThreadPostSerializer, ThreadCommentSerializer, ThreadSerializer, \
    ThreadUserSerializer, ThreadNestedSerializer, ThreadPostNestedSerializer, \
    ThreadFilterSerializer, ThreadUserInviteSerializer, \
    ThreadUserInviteUpdateSerializer, ThreadUserUpdateSerializer
from ..filters import ThreadFilter, MyThreadFilter

User = get_user_model()


class ThreadViewSet(MultiSerializerMixin,
                    viewsets.ReadOnlyModelViewSet):

    """Threads

    ### Routes ###

    * **[GET] /threads/**: List of threads
    * **[GET] /threads/<id>**: Detail of thread

    ### Additional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **doi=(str)**: Fetch thread associated with DOI
    * **pmi=(str)**: Fetch thread associated with PubMed id
    * **arx=(str)**: Fetch thread associated with Arxiv
    * **pii=(str)**: Fetch thread associated with PII Elsevier id
    * **type=(int)**: Fetch thread associated with type of thread
    * **title=(str)**: Fetch thread on title
    * **min_date=(str)**: Fetch thread published after min_date
    * **max_date=(str)**: Fetch thread published before max_date
    * **time_span=(str)**: Fetch thread published during last time_span days
    * **search=(str)**: Fetch thread on title, content, author first and last names
    * **third_party=(str)**: Fetch threads from third party (e.g 'pubpeer') or from etalia only ('false') (default: None)

    """

    queryset = Thread.objects.filter(privacy=THREAD_PUBLIC)
    serializer_class = {
        'default': ThreadSerializer,
        'nested': ThreadNestedSerializer,
    }
    permission_classes = (permissions.AllowAny,
                          )
    filter_backends = (DisabledHTMLContextualFilterBackend,
                       DisabledHTMLSearchFilterBackend,
                       )
    filter_class = ThreadFilter
    search_fields = ('title',
                     'content',
                     'paper__id_doi',
                     'user__first_name',
                     'user__last_name')


class MyThreadViewSet(ThreadViewSet,
                      viewsets.ModelViewSet):

    """
    My Threads

    ### Routes ###

    * **[GET, POST] /threads/**: List of threads
    * **[GET] /threads/filters**: Filter list for request threads list
    * **[GET, PUT, PATCH] /threads/<id\>/**: Thread instance
    * **[PUT, PATCH] /threads/<id\>/publish**: Publish Thread
    * **[GET] /threads/<id\>/neighbors**: Thread neighbors

    ### Additional Kwargs ###

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
    * **time_span=(int)**: Fetch only threads published in the past time_span days
    * **sort_by=(str)**: Sort threads by: 'date', 'score', 'published-date' (default = published-date)
    * **user_id[]=(int)**: Filter threads by user_id
    * **type[]=(int)**: Filter threads by type
    * **search=(str)**: Filter threads on title, owner first and last names
    * **doi=(str)**: Filter threads based on paper DOI

    ** Sub-routes: **

    * **[GET] /threads/<id\>/neighbors**: Thread neighbors:

        * **time_span=(int)**: Fetch neighbors threads published in the past time_span days (default=60)
    """

    queryset = Thread.objects.all()
    serializer_class = {
        'default': ThreadSerializer,
        'nested': ThreadNestedSerializer,
    }
    exclude_action_serializers = {
        'list': ['nested'],
    }
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
                          IsOwnerOrReadOnly,
                          ThreadIsNotYetPublishedIsOwnerIfDeleteMethod)

    filter_backends = (DisabledHTMLContextualFilterBackend,
                       DisabledHTMLSearchFilterBackend,
                       )
    filter_class = MyThreadFilter
    search_fields = ('title',
                     'content',
                     'paper__id_doi',
                     'user__first_name',
                     'user__last_name')

    SIZE_MAX_USER_FILTER = 40
    NEIGHBORS_TIME_SPAN = 60
    TIME_SPAN_DEFAULT = 30

    def get_thread_id(self):
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
            self.request.session['feeds-control-states'] = {
                'time_span': time_span,
                'search': search,
                'pin': 1 if pin == '1' else 0
            }
        else:
            self.request.session['threads-control-states'] = {
                'time_span': None,
                'search': search,
                'pin': 0
            }

    def get_queryset(self):
        # Store session control
        self.store_controls(self.request)
        # Defined baseline queryset
        self.queryset = Thread.objects.filter(
            ~(Q(published_at=None) & ~Q(user=self.request.user)),
            ~(Q(privacy=THREAD_PRIVATE) &
              ~(Q(threaduser__user=self.request.user) &
                Q(threaduser__participate=THREAD_JOINED))
              )
        )
        queryset = super(MyThreadViewSet, self).get_queryset()
        return self.get_serializer_class()\
            .setup_eager_loading(queryset, user=self.request.user)

    @detail_route(methods=['put', 'patch'],
                  permission_classes=(ThreadIsNotYetPublished,))
    def publish(self, request, pk=None):
        instance = get_object_or_404(self.queryset, pk=pk)
        self.check_object_permissions(request, instance)
        instance.publish()
        instance.embed()
        serializer = self.get_serializer(instance)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @detail_route(methods=['get'],
                  permission_classes=(ThreadIsPublished, ))
    def neighbors(self, request, pk=None):
        time_span = int(self.request.query_params.get(
            'time_span',
            settings.THREADS_DEFAULT_NEIGHBORS_TIMESPAN))
        instance = get_object_or_404(self.queryset, pk=pk)
        self.check_object_permissions(request, instance)
        neighbors = instance.get_neighbors(time_span)
        serializer = self.get_serializer(neighbors, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

    @list_route(methods=['get'])
    def filters(self, request):
        queryset = super(ThreadViewSet, self).get_queryset()
        queryset = queryset.select_related('user')
        data = list(self.filter_queryset(queryset))
        du = [d.user for d in data]

        us_count = Counter(du).most_common()
        users = []
        for u, c in us_count[:self.SIZE_MAX_USER_FILTER]:
            u.count = c
            if not u.first_name and not u.last_name:
                u.last_name = 'Anonymous User'
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

    ### Additional Kwargs ###

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
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
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
            if thread_id:
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

    ### Additional Kwargs ###

    ** List: **

    * **post_id=(int)**: Filter comments related to Post post_id

    """

    queryset = ThreadComment.objects.filter()
    serializer_class = {
        'default': ThreadCommentSerializer,
    }
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
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
            if post_id:
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
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
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

    ### Additional Kwargs ###

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
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
                          IsToUserOrOwnersReadOnly,
                          )

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            queryset = ThreadUserInvite.objects.filter(
                Q(from_user=self.request.user) |
                Q(to_user=self.request.user))
            # filter from_user
            from_user = self.request.query_params.get('from-user', None)
            if from_user:
                queryset = queryset.filter(from_user=from_user)
            # filter to_user
            to_user = self.request.query_params.get('to-user', None)
            if to_user:
                queryset = queryset.filter(to_user=to_user)
            # filter status
            status_ = self.request.query_params.get('status', None)
            if status_:
                queryset = queryset.filter(status=status_)

            return queryset

        return ThreadUserInvite.objects.all()

    def perform_create(self, serializer):
        serializer.save(from_user=self.request.user)