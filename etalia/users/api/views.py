# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db.models import Q

from rest_framework.decorators import detail_route, list_route
from rest_framework.response import Response
from rest_framework import viewsets, permissions, mixins, status
from rest_framework.exceptions import ParseError

from etalia.core.api.permissions import IsReadOnlyRequest, IsOwnerOrReadOnly, \
    IsInRelationship

from etalia.core.api.permissions import IsOwner, IsOwnerIfViewFull
from etalia.core.api.mixins import MultiSerializerMixin, One2OneNestedLinkSwitchMixin
from etalia.library.api.serializers import PaperSerializer, \
    PaperNestedSerializer

from etalia.core.api.mixins import MultiSerializerMixin
from .serializers import UserLibSerializer, UserSerializer, UserFullSerializer, \
    UserLibNestedSerializer, UserLibPaperSerializer, RelationshipSerializer

from ..models import UserLib, User, Relationship, UserLibPaper
from etalia.users.constants import RELATIONSHIP_BLOCKED, RELATIONSHIP_FOLLOWING, \
    RELATIONSHIP_STATUSES

class UserViewSet(MultiSerializerMixin,
                  mixins.ListModelMixin,
                  mixins.RetrieveModelMixin,
                  mixins.UpdateModelMixin,
                  viewsets.GenericViewSet):
    """
    Users

    ### Routes ###

    * **[GET] /users/**: List of user
    * **[GET] /users/<id\>/**: User instance
    * **[GET] /users/<id\>/followers/**: List of followers users
    * **[GET] /users/<id\>/following/**: List of following users
    * **[GET] /users/<id\>/blocked/**: List of blocked users

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'full',

    """

    queryset = User.objects.all()
    permission_classes = (permissions.IsAuthenticated,
                          IsOwnerOrReadOnly,
                          IsOwnerIfViewFull)
    serializer_class = {
        'default': UserSerializer,
        'full': UserFullSerializer,
    }
    exclude_action_serializers = {
        'list': 'full',
    }

    def render_list(self, queryset):
        page = self.paginate_queryset(queryset)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)

    @detail_route(permission_classes=(IsOwner, ))
    def following(self, request, pk=None):
        return self.render_list(self.request.user.following)

    @detail_route(permission_classes=(IsOwner, ))
    def followers(self, request, pk=None):
        return self.render_list(self.request.user.followers)

    @detail_route(permission_classes=(IsOwner, ))
    def blocked(self, request, pk=None):
        return self.render_list(self.request.user.blocked)


class UserLibViewSet(MultiSerializerMixin,
                     viewsets.ReadOnlyModelViewSet):
    """

    User library (currently only 1 per user)

    ### Routes ###

    * **[GET] /user-libs/**: List of papers in user libraries
    * **[GET] /user-libs/<id\>/**: User Library instance
    * **[GET] /user-libs/<id\>/papers**: List of papers in user library

    ### Optional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',

    """

    queryset = UserLib.objects.all()
    serializer_class = {
        'default': UserLibSerializer,
        'nested': UserLibNestedSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsReadOnlyRequest,
                          IsOwner)

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return UserLib.objects.filter(user=self.request.user)
        return UserLib.objects.all()

    @detail_route(methods=['GET'])
    def papers(self, request, pk=None):

        lib = self.get_object()
        page = self.paginate_queryset(lib.papers.all())
        self.serializer_class = {
            'default': PaperSerializer,
            'nested': PaperNestedSerializer,
        }

        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        serializer = self.get_serializer(lib.papers.all(), many=True)
        return Response(serializer.data)


class UserLibPaperViewSet(MultiSerializerMixin,
                          viewsets.ReadOnlyModelViewSet):
    """
    User Lib Paper: Relational Table between a [User Lib][ref1] and a [Paper][ref2]

    ### Routes ###

    * **[GET] /user-lib-papers/**: List of User Lib Paper
    * **[GET] /user-lib-papers/<id\>/**: User Lib Paper instance

    ### Optional Kwargs ###

    ** Detail: **

    * **view=(str)**: Reformat output. choices: 'nested',

    ** List: **

    * **view=(str)**: Reformat output. choices: 'nested',

    [ref1]: /api/v1/user/user-libs/
    [ref2]: /api/v1/library/papers/

    """

    queryset = UserLibPaper.objects.all()
    serializer_class = {
        'default': UserLibPaperSerializer,
        'nested': UserLibNestedSerializer
    }
    permission_classes = (permissions.IsAuthenticated,
                          IsReadOnlyRequest,
                          IsOwner)

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return UserLibPaper.objects.filter(userlib=self.request.user.lib)
        return UserLibPaper.objects.all()


class RelationshipViewSet(viewsets.ModelViewSet):
    """
    Relationship: between 2 users

    ### Routes ###

    * **[GET] /relationships/**: List of relationship for user
    * **[GET] /relationships/<id\>/**: Relationship instance for user

    ### Optional Kwargs ###

    ** List: **

    * **status=(int)**: Fetch only corresponding status
    * **from-user=(int)**: Fetch only relationships from user <id>
    * **to-user=(int)**: Fetch only relationships to user <id>

    """
    queryset = Relationship.objects.all()
    serializer_class = RelationshipSerializer
    permissions_classes = (permissions.IsAuthenticated,
                           IsInRelationship)

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

        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            queryset = Relationship.objects.\
                filter(Q(from_user=self.request.user) |
                       Q(to_user=self.request.user))

            if 'status' in self.request.query_params.keys():
                param = self.request.query_params.get('status', None)
                queryset = queryset.filter(status=param)
            if 'from_user' in self.request.query_params.keys():
                param = self.request.query_params.get('from-user', None)
                queryset = queryset.filter(from_user=param)
            if 'to_user' in self.request.query_params.keys():
                param = self.request.query_params.get('to-user', None)
                queryset = queryset.filter(to_user=param)

            return queryset
        return Relationship.objects.all()

    def perform_create(self, serializer):
        serializer.save(from_user=self.request.user)