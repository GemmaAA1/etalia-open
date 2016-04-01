# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import get_object_or_404

from rest_framework import viewsets, permissions, mixins, status
from rest_framework.response import Response
from rest_framework.decorators import detail_route
from rest_condition import And, Or, Not

from etalia.core.api.permissions import IsReadOnlyRequest, IsPostRequest, \
    IsDeleteRequest, IsPutPatchRequest, IsThreadMember, IsOwner, \
    IsJoinAction, IsLeaveAction, IsPinBanAction, IsStateAction, \
    IsNOTThreadMember, IsNOTStateAction

from ..models import Thread, ThreadPost, ThreadComment, ThreadUser


from .serializers import \
    ThreadPostSerializer, ThreadCommentSerializer, ThreadSerializer, \
    ThreadUserSerializer, ThreadNestedSerializer, ThreadPostNestedSerializer, \
    ThreadCommentNestedSerializer
from .mixins import ListRetrieveNestedMixin


class ThreadViewSet(ListRetrieveNestedMixin,
                    mixins.CreateModelMixin,
                    mixins.ListModelMixin,
                    mixins.RetrieveModelMixin,
                    mixins.UpdateModelMixin,
                    viewsets.GenericViewSet):
    """Thread view set

    Destroy (DELETE) routes is not provided
    """

    queryset = Thread.objects.all()
    serializer_class = ThreadSerializer
    serializer_nested_class = ThreadNestedSerializer
    permission_classes = (And(permissions.IsAuthenticated,
                              Or(And(IsReadOnlyRequest, ),
                                 And(IsPostRequest, ),
                                 And(IsPutPatchRequest, IsOwner, IsNOTStateAction),
                                 And(IsPutPatchRequest, IsStateAction),
                                 ),
                              ),
                          )

    def get_thread_id(self):
        return self.kwargs['pk']

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def perform_threaduser_action(self, request, action):
        """Perform ThreadUser related action (join, pin, etc.)"""
        thread = self.get_object()
        self.check_object_permissions(request, thread)
        instance, new = ThreadUser.objects.get_or_create(
            user=request.user,
            thread_id=self.get_thread_id())
        method = getattr(instance, action)
        method()
        # Return
        if new:
            return Response({}, status=status.HTTP_201_CREATED)
        else:
            return Response({}, status=status.HTTP_200_CREATED)

    @detail_route(methods=['patch', 'post'],
                  permission_classes=(IsNOTThreadMember, ))
    def join(self, request, pk=None):
        return self.perform_threaduser_action(request, 'join')

    @detail_route(methods=['patch', 'post'],
                  permission_classes=(IsThreadMember, ))
    def leave(self, request, pk=None):
        return self.perform_threaduser_action(request, 'leave')

    @detail_route(methods=['patch', 'post'])
    def pin(self, request, pk=None):
        return self.perform_threaduser_action(request, 'pin')

    @detail_route(methods=['patch', 'post'])
    def ban(self, request, pk=None):
        return self.perform_threaduser_action(request, 'ban')


class ThreadPostViewSet(ListRetrieveNestedMixin, viewsets.ModelViewSet):

    queryset = ThreadPost.objects.all()
    serializer_class = ThreadPostSerializer
    serializer_nested_class = ThreadPostNestedSerializer
    permission_classes = (And(permissions.IsAuthenticated,
                              Or(And(IsReadOnlyRequest, IsThreadMember),
                                 And(IsPostRequest, IsThreadMember),
                                 And(IsPutPatchRequest, IsThreadMember, IsOwner))), )

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def get_queryset(self):
        threads_joined = ThreadUser.objects\
            .filter(user=self.request.user, is_joined=True)\
            .values('thread')
        queryset = ThreadPost.objects.filter(thread__in=threads_joined)

        # filter based on ?thread_id
        thread_id = self.request.query_params.get('thread_id', None)
        if thread_id is not None:
            queryset = queryset.filter(thread_id=thread_id)

        return queryset


class ThreadCommentViewSet(ListRetrieveNestedMixin, viewsets.ModelViewSet):

    queryset = ThreadComment.objects.filter()
    serializer_class = ThreadCommentSerializer
    serializer_nested_class = ThreadCommentNestedSerializer
    permission_classes = (And(permissions.IsAuthenticated,
                              Or(And(IsReadOnlyRequest, IsThreadMember),
                                 And(IsPostRequest, IsThreadMember),
                                 And(IsPutPatchRequest, IsThreadMember, IsOwner))), )

    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

    def get_queryset(self):
        threads_joined = ThreadUser.objects\
            .filter(user=self.request.user, is_joined=True)\
            .values('thread')
        return ThreadComment.objects.filter(post__thread__in=threads_joined)


class ThreadUserViewSet(viewsets.ModelViewSet):

    queryset = ThreadUser.objects.filter()
    serializer_class = ThreadUserSerializer
    permission_classes = (And(permissions.IsAuthenticated, IsReadOnlyRequest), )

    def get_queryset(self):
        return ThreadUser.objects.filter(user=self.request.user)

