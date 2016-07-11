# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib.auth import get_user_model

from rest_framework import viewsets, permissions, mixins

from etalia.core.api.permissions import IsOwner

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
    permission_classes = (permissions.IsAuthenticated,
                          IsOwner,
                          )

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return UserPopOver.objects.filter(user=self.request.user)
        return UserPopOver.objects.all()


class PopOverViewSet(viewsets.ModelViewSet):
    """
    PopOver

    ### Routes ###

    * **[GET, POST] /popovers/**: List of PopOvers
    * **[GET, PUT, PATCH, DELETE] /popovers/<id\>/**: PopOver detail

    """

    queryset = PopOver.objects.all()
    serializer_class = PopOverSerializer
    permission_classes = (permissions.IsAuthenticated, )

    def perform_create(self, serializer):
        instance = serializer.save()
        # Bulk create userpopover
        us = User.objects.all()
        for user in us:
            UserPopOver.objects.create(popover=instance, user=user)

