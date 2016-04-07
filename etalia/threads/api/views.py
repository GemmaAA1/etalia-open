# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import viewsets, permissions, mixins, status
from rest_framework.response import Response
from rest_framework.decorators import detail_route
from rest_condition import And, Or

from etalia.core.api.permissions import IsThreadMember, IsOwner, \
    IsOwnerOrReadOnly, IsNOTThreadMember
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser
from ..constant import THREAD_JOINED, THREAD_LEFT, THREAD_PINNED, THREAD_BANNED
from .serializers import \
    ThreadPostSerializer, ThreadCommentSerializer, ThreadSerializer, \
    ThreadUserSerializer, ThreadNestedSerializer, ThreadPostNestedSerializer, \
    ThreadCommentNestedSerializer, PatchSerializer
from etalia.core.api.mixins import MultiSerializerMixin


class ThreadViewSet(MultiSerializerMixin,
                    mixins.CreateModelMixin,
                    mixins.ListModelMixin,
                    mixins.RetrieveModelMixin,
                    mixins.UpdateModelMixin,
                    viewsets.GenericViewSet):
    """
    Threads

    ### Routes ###

    * [GET, POST] /threads/: List of threads
    * [GET, PUT, PATCH] /threads/<id\>/: Thread instance

    ### Optional Kwargs ###

    ** Detail: **

    * ?view=(str): Reformat output. choices: 'nested',

    ** List: **

    * ?view=(str): Reformat output. choices: 'nested',
    * ?pinned=(int): Fetch only **pinned** threads for logged user if 1 (default = 0)
    * ?joined=(int): Fetch only **joined** threads for logged user if 1 (default = 0)
    * ?left=(int): Fetch only **left** threads for logged user if 1 (default = 0)
    * ?banned=(int): Fetch only **banned** threads for logged user if 1 (default = 0)

    """

    queryset = Thread.objects.all()
    serializer_class = {
        'default': ThreadSerializer,
        'nested': ThreadNestedSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsOwnerOrReadOnly)

    def get_thread_id(self):
        return self.kwargs['pk']

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def get_queryset(self):
        queryset = Thread.objects.all()

        # filter joined threads for user
        joined = self.request.query_params.get('joined', False)
        if joined:
            queryset = queryset.filter(threaduser__user=self.request.user,
                                       threaduser__participate=THREAD_JOINED)
        # filter pinned threads for user
        pinned = self.request.query_params.get('pinned', False)
        if pinned:
            queryset = queryset.filter(threaduser__user=self.request.user,
                                       threaduser__watch=THREAD_PINNED)
        # filter left threads for user
        left = self.request.query_params.get('left', False)
        if left:
            queryset = queryset.filter(threaduser__user=self.request.user,
                                       threaduser__participate=THREAD_LEFT)

        # filter banned threads for user
        banned = self.request.query_params.get('banned', False)
        if banned:
            queryset = queryset.filter(threaduser__user=self.request.user,
                                       threaduser__participate=THREAD_BANNED)
        return queryset


class ThreadPostViewSet(MultiSerializerMixin,
                        viewsets.ModelViewSet):

    """
    Thread Post: Post on thread

    ### Routes ###

    * [GET, POST] /posts/: List of posts
    * [GET, PUT, PATCH] /posts/<id\>/: Post instance

    ### Optional Kwargs ###

    ** Detail:**

    * ?view=(str): Reformat output. choices: 'nested',

    ** List: **

    * ?view=(str): Reformat output. choices: 'nested',
    * ?thread_id=(int): Filter comments related to Thread thread_id

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
            threads_joined = ThreadUser.objects\
                .filter(user=self.request.user, participate=True)\
                .values('thread')
            queryset = ThreadPost.objects.filter(thread__in=threads_joined)

            # filter based on ?thread_id
            thread_id = self.request.query_params.get('thread_id', None)
            if thread_id is not None:
                queryset = queryset.filter(thread_id=thread_id)

            return queryset
        return ThreadPost.objects.all()


class ThreadCommentViewSet(MultiSerializerMixin, viewsets.ModelViewSet):

    """
    Thread Comment: Comments on thread post

    ### Routes ###

    * [GET, POST] /comments/: List of comments
    * [GET, PUT, PATCH] /comments/<id\>/: Comment instance

    ### Optional Kwargs ###

    ** Detail:**

    * ?view=(str): Reformat output. choices: 'nested',

    ** List: **

    * ?view=(str): Reformat output. choices: 'nested',
    * ?post_id=(int): Filter comments related to Post post_id

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
            threads_joined = ThreadUser.objects\
                .filter(user=self.request.user, participate=THREAD_JOINED)\
                .values('thread')

            queryset = ThreadComment.objects.filter(post__thread__in=threads_joined)

            # filter based on ?post_id
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

    * [GET, POST] /states/: List of Thread User states
    * [GET, PUT, PATCH] /states/<id\>/: Thread User state instance

    **Deprecated**:

    * [PUT, PATCH] /states/<id\>/join: User **join** thread
    * [PUT, PATCH] /states/<id\>/leave: User **leave** thread
    * [PUT, PATCH] /states/<id\>/pin: User **pin** thread
    * [PUT, PATCH] /states/<id\>/ban: User **ban** thread

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
                  permission_classes=(IsNOTThreadMember, ))
    def join(self, request, pk=None):
        return self.perform_threaduser_action(request, 'join')

    @detail_route(methods=['put', 'patch'],
                  permission_classes=(IsThreadMember, ))
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



