# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from etalia.users.serializers import UserSerializer
from etalia.library.serializers import PaperSerializer

from .models import Thread, ThreadPost, ThreadPostComment
from etalia.users.models import UserThread


class UserThreadSerializer(serializers.ModelSerializer):

    class Meta:
        model = UserThread
        fields = ('id',
                  'is_pinned',
                  'is_banned',
                  'is_joined',
                  'is_left',
                  'first_joined_at',
                  'last_left_at',
                  'num_comments')


class UserThreadCountSerializer(serializers.ModelSerializer):

    joined_count = serializers.SerializerMethodField()
    pinned_count = serializers.SerializerMethodField()
    left_count = serializers.SerializerMethodField()

    class Meta:
        model = UserThread
        fields = ('id',
                  'is_pinned',
                  'is_banned',
                  'is_joined',
                  'is_left',
                  'first_joined_at',
                  'last_left_at',
                  'num_comments')

    def get_joined_count(self, obj):
        return obj.user.userthread_set.filter(is_joined=True).count()

    def get_left_count(self, obj):
        return obj.user.userthread_set.filter(is_left=True).count()

    def get_pinned_count(self, obj):
        return obj.user.userthread_set.filter(is_pinned=True).count()


class ThreadPostCommentSerializer(serializers.ModelSerializer):

    author = UserSerializer(many=False, read_only=True)

    class Meta:
        model = ThreadPostComment
        fields = ('id',
                  'content',
                  'author',
                  'created',
                  'modified',
                  'post')


class ThreadPostSerializer(serializers.ModelSerializer):

    author = UserSerializer(many=False, read_only=True)
    comments = ThreadPostCommentSerializer(many=True, read_only=True)

    class Meta:
        model = ThreadPost
        fields = ('id',
                  'content',
                  'author',
                  'created',
                  'modified',
                  'comments',
                  'thread')


class ThreadSerializer(serializers.ModelSerializer):

    members = UserSerializer(many=True, read_only=True)
    owner = UserSerializer(many=False, read_only=True)
    paper = PaperSerializer(many=False, read_only=True)
    posts = ThreadPostSerializer(many=True, read_only=True)
    state = UserThreadSerializer(many=False, read_only=True)

    class Meta:
        model = Thread
        fields = ('id',
                  'type',
                  'title',
                  'owner',
                  'privacy',
                  'paper',
                  'title',
                  'content',
                  'members',
                  'posts',
                  'created',
                  'modified')


class ThreadContextualizedSerializer(serializers.ModelSerializer):

    members = UserSerializer(many=True, read_only=True)
    owner = UserSerializer(many=False, read_only=True)
    paper = PaperSerializer(many=False, read_only=True)
    posts = ThreadPostSerializer(many=True, read_only=True)
    state = UserThreadSerializer(many=False, read_only=True)

    class Meta:
        model = Thread
        fields = ('id',
                  'type',
                  'title',
                  'owner',
                  'privacy',
                  'paper',
                  'title',
                  'content',
                  'members',
                  'posts',
                  'created',
                  'modified',
                  'state')

    def get_state(self, obj):
        return UserThread(user=self.context['request'].user,
                          thread=obj)



