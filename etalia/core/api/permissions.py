# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import permissions

from etalia.threads.models import Thread, ThreadPost, ThreadComment


class IsReadOnlyRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method in permissions.SAFE_METHODS


class IsPostRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method == "POST"


class IsPutPatchRequest(permissions.BasePermission):

    def has_permission(self, request, view):
        return request.method in ["PUT", "PATCH"]


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


class IsOwner(permissions.BasePermission):
    """
    Custom permission to only allow owners of an object to edit it
    GIVEN THAT user is an attribute of the object pointing to owner
    """

    def has_object_permission(self, request, view, obj):
        return obj.user == request.user
