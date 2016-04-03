# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import viewsets

from etalia.core.api.permissions import IsReadOnlyRequest
from etalia.core.api.mixins import MultiSerializerMixin

from ..models import Paper, Author, Journal

from .serializers import PaperSerializer, JournalSerializer, AuthorSerializer, \
    PaperNestedSerializer


class PaperViewSet(MultiSerializerMixin, viewsets.ReadOnlyModelViewSet):

    queryset = Paper.objects.all()
    serializer_class = {
        'default': PaperSerializer,
         'nested': PaperNestedSerializer
    }
    permission_classes = (IsReadOnlyRequest, )


class JournalViewSet(viewsets.ReadOnlyModelViewSet):

    queryset = Journal.objects.all()
    serializer_class = JournalSerializer
    permission_classes = (IsReadOnlyRequest, )


class AuthorViewSet(viewsets.ReadOnlyModelViewSet):

    queryset = Author.objects.all()
    serializer_class = AuthorSerializer
    permission_classes = (IsReadOnlyRequest, )