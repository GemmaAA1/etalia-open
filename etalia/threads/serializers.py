# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from etalia.users.serializers import UserSerializer
from etalia.library.serializers import PaperSerializer

from .models import Thread, ThreadPost, ThreadPostComment


class ThreadPostCommentSerializer(serializers.ModelSerializer):

    author = UserSerializer(many=False, read_only=True)

    class Meta:
        model = ThreadPostComment
        fields = ('content',
                  'author',
                  'created',
                  'modified',
                  'post')


class ThreadPostSerializer(serializers.ModelSerializer):

    author = UserSerializer(many=False, read_only=True)
    comments = ThreadPostCommentSerializer(many=True, read_only=True)

    class Meta:
        model = ThreadPost
        fields = ('content',
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

    class Meta:
        model = Thread
        fields = ('type',
                  'title',
                  'owner',
                  'privacy',
                  'members',
                  'paper',
                  'title',
                  'content',
                  'posts')




