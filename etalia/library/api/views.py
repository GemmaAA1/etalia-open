# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import viewsets

from etalia.core.api.permissions import IsReadOnlyRequest
from etalia.core.api.mixins import MultiSerializerMixin

from ..models import Paper, Author, Journal

from .serializers import PaperSerializer, JournalSerializer, AuthorSerializer, \
    PaperNestedSerializer


class PaperViewSet(MultiSerializerMixin, viewsets.ReadOnlyModelViewSet):
    """
    Paper

    ### Routes ###

    * [GET] /papers/: List of papers
    * [GET] /papers/<id\>/: Paper instance
    """

    queryset = Paper.objects.all()
    serializer_class = {
        'default': PaperSerializer,
         'nested': PaperNestedSerializer
    }
    permission_classes = (IsReadOnlyRequest, )


class JournalViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Paper

    ### Routes ###

    * [GET] /journals/: List of journals
    * [GET] /journals/<id\>/: Journal instance
    """
    queryset = Journal.objects.all()
    serializer_class = JournalSerializer
    permission_classes = (IsReadOnlyRequest, )


class AuthorViewSet(viewsets.ReadOnlyModelViewSet):
    """
    Author

    ### Routes ###

    * [GET] /authors/: List of authors
    * [GET] /authors/<id\>/: Author instance
    """
    queryset = Author.objects.all()
    serializer_class = AuthorSerializer
    permission_classes = (IsReadOnlyRequest, )