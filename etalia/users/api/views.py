# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from rest_framework import viewsets, permissions

from etalia.core.api.permissions import IsOwner

from .serializers import UserLibSerializer, UserSerializer

from ..models import UserLib, User


class UserViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Returns a list of users

    """

    queryset = User.objects.all()
    serializer_class = UserSerializer
    permissions_classes = ()


class UserLibViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Returns a list of user libraries (currenlty only per user)

    ### Optional Kwargs ###

    ** All: **

    * view=(str): Reformat output. choices: 'nested',

    """

    queryset = UserLib.objects.all()
    serializer_class = UserLibSerializer
    permission_classes = (IsOwner, )

    def get_queryset(self):
        return UserLib.objects.filter(user=self.request.user)
