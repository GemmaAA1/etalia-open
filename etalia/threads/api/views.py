# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import viewsets, permissions, mixins, status
from rest_framework.response import Response
from rest_framework.decorators import detail_route
from rest_condition import And, Or, Not

from etalia.core.api.permissions import IsReadOnlyRequest, IsPostRequest, \
    IsDeleteRequest, IsPutPatchRequest, IsThreadMember, IsOwner

from ..models import Thread, ThreadPost, ThreadComment, ThreadUser


from .serializers import FullThreadSerializer, BasicThreadSerializer, \
    ThreadPostSerializer, ThreadCommentSerializer, \
    CreateThreadSerializer, UpdateThreadSerializer, ThreadUserSerializer, \
    CreateUpdateThreadPostSerializer


class ThreadViewSet(mixins.CreateModelMixin,
                    mixins.ListModelMixin,
                    mixins.RetrieveModelMixin,
                    mixins.UpdateModelMixin,
                    viewsets.GenericViewSet):
    """Thread view set

    Destroy (DELETE) routes is not provided
    """

    queryset = Thread.objects.all()
    serializer_class = None   # action/permission dependent
    permission_classes = (And(permissions.IsAuthenticated,
                              Or(And(IsReadOnlyRequest, ),
                                 And(IsPostRequest, ),
                                 And(IsPutPatchRequest, IsOwner))), )

    def get_serializer_class(self):
        # Read only method
        if IsReadOnlyRequest().has_permission(self.request, self):
            if self.action == 'retrieve' \
                and IsThreadMember().has_object_permission(self.request,
                                                           self,
                                                           self.get_object()):
                return FullThreadSerializer
            return BasicThreadSerializer
        elif IsPostRequest().has_permission(self.request, self):
            return CreateThreadSerializer
        elif IsPutPatchRequest().has_permission(self.request, self):
            return UpdateThreadSerializer
        return super(ThreadViewSet, self).get_serializer_class()

    def get_thread_id(self):
        return self.kwargs['pk']

    def perform_action(self, request, action):
        """Perform ThreadUser related action (join, pin, etc.)"""
        instance, new = ThreadUser.objects.get_or_create(
            user=request.user,
            thread_id=self.get_thread_id())
        method = getattr(instance, action)
        method()
        return Response(ThreadUserSerializer(instance).data)

    @detail_route(methods=['patch', 'post'],
                  permission_classes=[Not(IsThreadMember)])
    def join(self, request, pk=None):
        return self.perform_action(request, 'join')

    @detail_route(methods=['patch', 'post'],
                  permission_classes=[IsThreadMember])
    def leave(self, request, pk=None):
        return self.perform_action(request, 'leave')

    @detail_route(methods=['patch', 'post'])
    def pin(self, request, pk=None):
        return self.perform_action(request, 'pin')

    @detail_route(methods=['patch', 'post'])
    def ban(self, request, pk=None):
        return self.perform_action(request, 'ban')


class ThreadPostViewSet(viewsets.ModelViewSet):

    queryset = ThreadPost.objects.all()
    serializer_class = ThreadPostSerializer
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
        return ThreadPost.objects.filter(thread__in=threads_joined)

    def update(self, request, *args, **kwargs):
        return super(ThreadPostViewSet, self).update(request, *args, **kwargs)

    # def get_serializer_class(self):
    #     if IsPostRequest().has_permission(self.request, self):
    #         return CreateUpdateThreadPostSerializer
    #     if IsPutPatchRequest().has_permission(self.request, self):
    #         return CreateUpdateThreadPostSerializer
    #     return super(ThreadPostViewSet, self).get_serializer_class()
    #
    # def create(self, request, *args, **kwargs):
    #     """Simplified serializer for user creation"""
    #     kwargs['context'] = self.get_serializer_context()
    #     serializer_in = CreateUpdateThreadPostSerializer(data=request.data,
    #                                                   **kwargs)
    #     serializer_in.is_valid(raise_exception=True)
    #     instance = serializer_in.save()
    #     serializer_out = self.get_serializer(instance=instance)
    #     headers = self.get_success_headers(serializer_out.data)
    #     return Response(serializer_out.data, status=status.HTTP_201_CREATED, headers=headers)
    
    # def update(self, request, *args, **kwargs):
    #     return super(ThreadPostViewSet, self).update(request, *args, **kwargs)


class ThreadCommentViewSet(viewsets.ModelViewSet):

    queryset = ThreadComment.objects.filter()
    serializer_class = ThreadCommentSerializer
    permission_classes = (And(permissions.IsAuthenticated,
                              Or(And(IsReadOnlyRequest, IsThreadMember),
                                 And(IsPostRequest, IsThreadMember),
                                 And(IsPutPatchRequest, IsThreadMember, IsOwner))), )

    def get_queryset(self):
        threads_joined = ThreadUser.objects\
            .filter(user=self.request.user, is_joined=True)\
            .values('thread')
        return ThreadComment.objects.filter(post__thread__in=threads_joined)

    def get_serializer_class(self):
        if IsPostRequest().has_permission(self.request, self):
            return CreateUpdateThreadCommentSerializer
        if IsPutPatchRequest().has_permission(self.request, self):
            return CreateUpdateThreadCommentSerializer
        return super(ThreadCommentViewSet, self).get_serializer_class()

