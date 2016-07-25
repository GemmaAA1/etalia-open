# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib.auth import get_user_model
from django.views.decorators.cache import never_cache

from rest_framework import viewsets, permissions, mixins
from rest_framework import filters

from etalia.core.api.permissions import IsOwner, IsAdminOrReadOnly, IsSessionAuthenticatedOrReadOnly

from .serializers import UserPopOverSerializer, PopOverSerializer
from ..models import UserPopOver, PopOver


User = get_user_model()


class PopOverStateViewSet(mixins.CreateModelMixin,
                          mixins.ListModelMixin,
                          mixins.UpdateModelMixin,
                          mixins.RetrieveModelMixin,
                          viewsets.GenericViewSet):
    """
    PopOver User state: Relation between a User and a PopOver

    ### Routes ###

    * **[GET, POST] /states/**: List of Paper/User states
    * **[GET, PUT, PATCH] /states/<id\>/**: Paper/User state instance

    """

    queryset = UserPopOver.objects.all()
    serializer_class = UserPopOverSerializer
    permission_classes = (IsSessionAuthenticatedOrReadOnly,
                          IsOwner,
                          )

    filter_backends = (filters.DjangoFilterBackend, )
    filter_fields = ('status', 'display')

    @never_cache
    def list(self, request, *args, **kwargs):
        return super(PopOverStateViewSet, self).list(request, *args, **kwargs)

    @never_cache
    def retrieve(self, request, *args, **kwargs):
        return super(PopOverStateViewSet, self).list(request, *args, **kwargs)

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return UserPopOver.objects\
                .filter(user=self.request.user)\
                .order_by('-popover__type', 'popover__priority')
        return UserPopOver.objects.all()


class PopOverViewSet(viewsets.ModelViewSet):
    """
    PopOver

    ### Routes ###

    * **[GET, POST] /popovers/**: List of PopOvers
    * **[GET, PUT, PATCH, DELETE] /popovers/<id\>/**: PopOver detail

    """

    queryset = PopOver.objects.all().order_by('-type', '-priority')
    serializer_class = PopOverSerializer
    permission_classes = (IsAdminOrReadOnly, )

    filter_backends = (filters.DjangoFilterBackend, )
    filter_fields = ('type', )

    def perform_create(self, serializer):
        instance = serializer.save()
        instance.init()
