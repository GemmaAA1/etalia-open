# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin
from etalia.library.api.serializers import PaperNestedSerializer
from etalia.threads.api.serializers import ThreadNestedSerializer
from etalia.users.api.serializers import UserSerializer

from ..models import StreamPapers, TrendPapers, ThreadFeedThreads, Stream, \
    Trend, ThreadFeed


class StreamSerializer(One2OneNestedLinkSwitchMixin,
                       serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Stream
        extra_kwargs = {
            'link': {'view_name': 'api:stream-detail'},
            'user': {'view_name': 'api:user-detail'},
            'papers': {'view_name': 'api:paper-detail'}
        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer},
        }


class StreamPapersSerializer(One2OneNestedLinkSwitchMixin,
                             serializers.HyperlinkedModelSerializer):

    class Meta:
        model = StreamPapers
        extra_kwargs = {
            'link': {'view_name': 'api:streampapers-detail'},
            'paper': {'view_name': 'api:paper-detail'},
            'stream': {'view_name': 'api:stream-detail'}
        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'paper': {'serializer': PaperNestedSerializer},
        }


class TrendSerializer(One2OneNestedLinkSwitchMixin,
                      serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Trend
        extra_kwargs = {
            'link': {'view_name': 'api:trend-detail'},
            'user': {'view_name': 'api:user-detail'},
            'papers': {'view_name': 'api:paper-detail'},
        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer},
        }


class TrendPapersSerializer(One2OneNestedLinkSwitchMixin,
                            serializers.HyperlinkedModelSerializer):

    class Meta:
        model = TrendPapers
        extra_kwargs = {
            'link': {'view_name': 'api:trendpapers-detail'},
            'paper': {'view_name': 'api:paper-detail'},
            'trend': {'view_name': 'api:trend-detail'}

        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'paper': {'serializer': PaperNestedSerializer},
        }


class ThreadFeedSerializer(One2OneNestedLinkSwitchMixin,
                           serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadFeed
        extra_kwargs = {
            'link': {'view_name': 'api:threadfeed-detail'},
            'user': {'view_name': 'api:user-detail'},
            'threads': {'view_name': 'api:thread-detail'}
        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer},
        }


class ThreadFeedThreadsSerializer(One2OneNestedLinkSwitchMixin,
                                  serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadFeedThreads
        extra_kwargs = {
            'link': {'view_name': 'api:threadfeedthreads-detail'},
            'thread': {'view_name': 'api:thread-detail'},
            'threadfeed': {'view_name': 'api:threadfeed-detail'}
        }
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'thread': {'serializer': ThreadNestedSerializer},
        }