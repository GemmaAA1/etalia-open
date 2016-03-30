# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from etalia.users.api.serializers import UserSerializer
from etalia.library.api.serializers import PaperSerializer
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser
from ..constant import THREAD_PRIVACIES, THREAD_TYPES
from .utils import KeyValueField


class ThreadUserSerializer(serializers.ModelSerializer):
    class Meta:
        model = ThreadUser
        fields = ('id',
                  'is_pinned',
                  'is_banned',
                  'is_joined',
                  'is_left',
                  'first_joined_at',
                  'last_left_at',
                  'num_comments')


class ThreadUserCountSerializer(serializers.ModelSerializer):
    joined_count = serializers.SerializerMethodField()
    pinned_count = serializers.SerializerMethodField()
    left_count = serializers.SerializerMethodField()

    class Meta:
        model = ThreadUser
        fields = ('id',
                  'is_pinned',
                  'is_banned',
                  'is_joined',
                  'is_left',
                  'first_joined_at',
                  'last_left_at',
                  'num_comments',
                  'joined_count',
                  'pinned_count',
                  'left_count')

    def get_joined_count(self, obj):
        return obj.user.ThreadUser_set.filter(is_joined=True).count()

    def get_left_count(self, obj):
        return obj.user.ThreadUser_set.filter(is_left=True).count()

    def get_pinned_count(self, obj):
        return obj.user.ThreadUser_set.filter(is_pinned=True).count()


class CreateUpdateThreadCommentSerializer(serializers.ModelSerializer):
    class Meta:
        model = ThreadComment
        fields = ('id',
                  'content',
                  'post')

    def save(self, **kwargs):
        return super(CreateUpdateThreadCommentSerializer, self).save(
            user=self.context['request'].user,
            **kwargs)


class ThreadCommentSerializer(serializers.ModelSerializer):
    user = UserSerializer(many=False, read_only=True)

    class Meta:
        model = ThreadComment
        fields = ('id',
                  'content',
                  'user',
                  'created',
                  'modified',
                  'post')


class CreateUpdateThreadPostSerializer(serializers.ModelSerializer):
    class Meta:
        model = ThreadPost
        fields = ('id',
                  'thread',
                  'content',)

    def create(self, validated_data):
        return super(CreateUpdateThreadPostSerializer, self).create(
            validated_data)

    def save(self):
        return super(CreateUpdateThreadPostSerializer, self).save(
            user=self.context['request'].user)


class ThreadPostSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadPost
        extra_kwargs = {'link': {'view_name': 'api:threadpost-detail'},
                        'comments': {'view_name': 'api:threadcomment-detail'},
                        'user': {'view_name': 'api:user-detail'},
                        'thread': {'view_name': 'api:thread-detail'}
                        }
        fields = (
            'link',
            'thread',
            'user',
            'content',
            'created',
            'modified',
            'comments')
        read_only_fields = (
            'link',
            'user',
            'created',
            'modified',
            'comments')


class CreateThreadSerializer(serializers.ModelSerializer):
    """Update serializer for Thread"""

    type = KeyValueField(labels=dict(THREAD_TYPES))
    privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))

    class Meta:
        model = Thread
        fields = ('id',
                  'type',
                  'paper',
                  'privacy',
                  'title')

    def save(self, **kwargs):
        return super(CreateThreadSerializer, self).save(
            user=self.context['request'].user,
            **kwargs)


class UpdateThreadSerializer(serializers.ModelSerializer):
    """Update serializer for Thread"""

    class Meta:
        model = Thread
        fields = ('id',
                  'title',
                  'content')


class BasicThreadSerializer(serializers.ModelSerializer):
    """Basic serializer for Thread when request.user is not a member"""
    user = UserSerializer(many=False, read_only=True)
    paper = PaperSerializer(many=False, read_only=True)
    state = serializers.SerializerMethodField(read_only=True)
    link = serializers.HyperlinkedIdentityField(view_name='api:thread-detail')
    type = KeyValueField(labels=dict(THREAD_TYPES))
    privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))

    class Meta:
        model = Thread
        fields = ('id',
                  'link',
                  'type',
                  'title',
                  'user',
                  'privacy',
                  'state',
                  'paper',
                  'title',
                  'content',
                  'created',
                  'modified')

    def get_state(self, obj):
        if ThreadUser.objects.filter(user=self.context['request'].user,
                                     thread=obj).exists():
            instance = ThreadUser.objects.get(user=self.context['request'].user,
                                              thread=obj)
            return ThreadUserSerializer(instance=instance,
                                        context={'request': self.context[
                                            'request']}).data
        else:
            return None


class FullThreadSerializer(serializers.ModelSerializer):
    """Full thread serializer when request.user is a member"""
    members = UserSerializer(many=True, read_only=True)
    user = UserSerializer(many=False, read_only=True)
    paper = PaperSerializer(many=False, read_only=True)
    posts = ThreadPostSerializer(many=True, read_only=True)
    state = serializers.SerializerMethodField(read_only=True)
    type = KeyValueField(labels=dict(THREAD_TYPES))
    privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))
    link = serializers.HyperlinkedIdentityField(view_name='api:thread-detail')

    class Meta:
        model = Thread
        fields = ('id',
                  'link',
                  'type',
                  'title',
                  'user',
                  'privacy',
                  'state',
                  'paper',
                  'title',
                  'content',
                  'members',
                  'posts',
                  'created',
                  'modified')

    def get_state(self, obj):
        if ThreadUser.objects.filter(user=self.context['request'].user,
                                     thread=obj).exists():
            instance = ThreadUser.objects.get(user=self.context['request'].user,
                                              thread=obj)
            return ThreadUserSerializer(instance=instance,
                                        context={'request': self.context[
                                            'request']}).data
        else:
            return None
