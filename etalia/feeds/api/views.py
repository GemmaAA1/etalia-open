# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.decorators.cache import never_cache

from rest_framework import viewsets, permissions
from rest_framework import filters

from .serializers import StreamPapersSerializer, TrendPapersSerializer, \
    ThreadFeedThreadsSerializer, StreamSerializer, TrendSerializer, \
    ThreadFeedSerializer

from ..models import StreamPapers, TrendPapers, ThreadFeedThreads, Stream, \
    Trend, ThreadFeed


class StreamViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = StreamSerializer
    queryset = Stream.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return Stream.objects\
                .filter(user=self.request.user)
        return Stream.objects.all()


class StreamPapersViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = StreamPapersSerializer
    queryset = StreamPapers.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return StreamPapers.objects\
                .filter(stream=self.request.user.streams.first())
        return StreamPapers.objects.all()


class TrendViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = TrendSerializer
    queryset = Trend.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return Trend.objects\
                .filter(user=self.request.user)
        return Trend.objects.all()


class TrendPapersViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = TrendPapersSerializer
    queryset = TrendPapers.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return TrendPapers.objects\
                .filter(trend=self.request.user.trends.first())
        return StreamPapers.objects.all()


class ThreadFeedViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = ThreadFeedSerializer
    queryset = ThreadFeed.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return ThreadFeed.objects\
                .filter(user=self.request.user)
        return ThreadFeed.objects.all()


class ThreadFeedThreadsViewSets(viewsets.ReadOnlyModelViewSet):

    serializer_class = ThreadFeedThreadsSerializer
    queryset = ThreadFeedThreads.objects.all()

    def get_queryset(self):
        # to raise proper 403 status code on not allowed access
        if self.action == 'list':
            return ThreadFeedThreads.objects\
                .filter(threadfeed=self.request.user.threadfeeds.first())
        return StreamPapers.objects.all()
