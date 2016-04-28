# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from django.contrib.auth import get_user_model
from rest_framework import permissions

from etalia.threads.models import Thread, ThreadPost, ThreadComment, \
    ThreadUserInvite
from etalia.users.models import UserLibPaper, Relationship

User = get_user_model()


class IsReadOnlyRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method in permissions.SAFE_METHODS


class IsPostRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "POST"


class IsPatchRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "PATCH"


class IsPutRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "PUT"


class IsPutRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "PUT"


class IsDeleteRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "DELETE"


class IsThreadMember(permissions.BasePermission):
    """
    Custom permission for user as member of thread
    """

    def has_object_permission(self, request, view, obj):
        if obj.__class__ == Thread:
            return request.user in obj.members
        if obj.__class__ == ThreadPost:
            return request.user in obj.thread.members
        if obj.__class__ == ThreadComment:
            return request.user in obj.post.thread.members


class IsNOTThreadMember(IsThreadMember):

    def has_object_permission(self, request, view, obj):
        return not super(IsNOTThreadMember, self).has_object_permission(request, view, obj)


class IsOwner(permissions.BasePermission):
    """
    Custom permission to only allow owners of an object to edit it
    GIVEN THAT user is an attribute of the object pointing to owner
    """

    def has_object_permission(self, request, view, obj):
        if obj.__class__ == UserLibPaper:
            return request.user == obj.userlib.user
        if obj.__class__ == User:
            return request.user == obj
        if obj.__class__ == ThreadUserInvite:
            return request.user in [obj.from_user, obj.to_user]
        if obj.__class__ == Relationship:
            return request.user == obj.from_user
        return obj.user == request.user


class IsOwnerOrReadOnly(IsOwner):
    """
    Custom permission to only allow owners of an object to edit it.
    """

    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        # Write permissions are only allowed to the owner of the snippet.
        return super(IsOwnerOrReadOnly, self).has_object_permission(request, view, obj)


class IsOwnerIfViewFull(IsOwner):
    """
    Custom permission to only allow owners of an object to edit it
    GIVEN THAT user is an attribute of the object pointing to owner
    """

    def has_object_permission(self, request, view, obj):
        if request.query_params.get('view', None) == 'full':
            return obj == request.user
        else:
            return True


class IsInRelationship(IsOwner):
    """
    Custom permission to only allow user involved in the relationship.
    """

    def has_object_permission(self, request, view, obj):
        if obj.__class__ == Relationship:
            return request.user in [obj.from_user, obj.to_user]


class IsStateAction(permissions.BasePermission):
    """
    Custom permission for ThreadUser action verbs
    """

    def has_permission(self, request, view):
        return view.action in ['join', 'leave', 'pin', 'ban']


class IsNOTStateAction(permissions.BasePermission):
    """
    Custom permission for ThreadUser action verbs
    """

    def has_permission(self, request, view):
        return not super(IsNOTStateAction, self).has_permission(request, view)


class IsJoinAction(permissions.BasePermission):
    """
    Custom permission for ThreadUser action verbs
    """

    def has_permission(self, request, view):
        return view.action == 'join'


class IsLeaveAction(permissions.BasePermission):
    """
    Custom permission for ThreadUser action verbs
    """

    def has_permission(self, request, view):
        return view.action == 'leave'


class IsPinBanAction(permissions.BasePermission):
    """
    Custom permission for ThreadUser action verbs
    """

    def has_permission(self, request, view):
        return view.action in ['pin', 'ban']


class ThreadIsNotYetPublished(permissions.BasePermission):
    """
    Custom permission to check if Thread has yet been published
    """

    def has_object_permission(self, request, view, obj):
        return obj.published_at is None


class ThreadIsNotYetPublishedIsOwnerIfDeleteMethod(IsOwner):
    """
    Custom permission to only allow owners of an object to edit it
    GIVEN THAT user is an attribute of the object pointing to owner
    """

    def has_object_permission(self, request, view, obj):
        if view.action in ['destroy']:
            return obj.published_at is None and obj == request.user
        return True


class ThreadIsPublished(permissions.BasePermission):
    """
    Custom permission to check if Thread has yet been published
    """

    def has_object_permission(self, request, view, obj):
        return obj.published_at is not None