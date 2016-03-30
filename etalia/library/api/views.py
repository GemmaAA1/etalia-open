# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import viewsets

from etalia.core.api.permissions import IsReadOnlyRequest

from ..models import Paper, Author, Journal

from .serializers import PaperSerializer, JournalSerializer


class PaperViewSet(viewsets.ReadOnlyModelViewSet):

    queryset = Paper.objects.all()
    serializer_class = PaperSerializer
    permission_classes = (IsReadOnlyRequest, )


class JournalViewSet(viewsets.ReadOnlyModelViewSet):

    queryset = Journal.objects.all()
    serializer_class = JournalSerializer
    permission_classes = (IsReadOnlyRequest, )
