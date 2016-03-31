# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf import settings

from rest_framework import serializers
from rest_framework.reverse import reverse

from etalia.users.api.serializers import UserSerializer
from etalia.library.api.serializers import PaperSerializer
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser
from ..constant import THREAD_PRIVACIES, THREAD_TYPES
from .utils import KeyValueField


class ThreadUserSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadUser
        extra_kwargs = {
            'link': {'view_name': 'api:threaduser-detail'},
        }
        fields = (
            'link',
            'is_pinned',
            'is_banned',
            'is_joined',
            'is_left',
            'first_joined_at',
            'last_left_at',
            'num_comments')
        read_only_fields = (
            'is_pinned',
            'is_banned',
            'is_joined',
            'is_left',
            'first_joined_at',
            'last_left_at',
            'num_comments'
        )


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


class ThreadCommentSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadComment
        extra_kwargs = {
            'link': {'view_name': 'api:threadcomment-detail'},
            'user': {'view_name': 'api:user-detail'},
            'post': {'view_name': 'api:threadpost-detail'}
        }
        fields = (
            'link',
            'user',
            'post',
            'content',
            'created',
            'modified',
        )
        read_only_fields = (
            'link',
            'user',
            'created',
            'modified')


class ThreadPostSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ThreadPost
        extra_kwargs = {
            'link': {'view_name': 'api:threadpost-detail'},
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


# class CreateThreadSerializer(serializers.ModelSerializer):
#     """Update serializer for Thread"""
#
#     type = KeyValueField(labels=dict(THREAD_TYPES))
#     privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))
#
#     class Meta:
#         model = Thread
#         fields = ('id',
#                   'type',
#                   'paper',
#                   'privacy',
#                   'title')
#
#     def save(self, **kwargs):
#         return super(CreateThreadSerializer, self).save(
#             user=self.context['request'].user,
#             **kwargs)
#
#
# class UpdateThreadSerializer(serializers.ModelSerializer):
#     """Update serializer for Thread"""
#
#     class Meta:
#         model = Thread
#         fields = ('id',
#                   'title',
#                   'content')


# class BasicThreadSerializer(serializers.ModelSerializer):
#     """Basic serializer for Thread when request.user is not a member"""
#     user = UserSerializer(many=False, read_only=True)
#     paper = PaperSerializer(many=False, read_only=True)
#     state = serializers.SerializerMethodField(read_only=True)
#     link = serializers.HyperlinkedIdentityField(view_name='api:thread-detail')
#     type = KeyValueField(labels=dict(THREAD_TYPES))
#     privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))
#
#     class Meta:
#         model = Thread
#         fields = ('id',
#                   'link',
#                   'type',
#                   'title',
#                   'user',
#                   'privacy',
#                   'state',
#                   'paper',
#                   'title',
#                   'content',
#                   'created',
#                   'modified')
#
#     def get_state(self, obj):
#         if ThreadUser.objects.filter(user=self.context['request'].user,
#                                      thread=obj).exists():
#             instance = ThreadUser.objects.get(user=self.context['request'].user,
#                                               thread=obj)
#             return ThreadUserSerializer(instance=instance,
#                                         context={'request': self.context[
#                                             'request']}).data
#         else:
#             return None


# class FullThreadSerializer(serializers.ModelSerializer):
#     """Full thread serializer when request.user is a member"""
#     members = UserSerializer(many=True, read_only=True)
#     user = UserSerializer(many=False, read_only=True)
#     paper = PaperSerializer(many=False, read_only=True)
#     posts = ThreadPostSerializer(many=True, read_only=True)
#     state = serializers.SerializerMethodField(read_only=True)
#     type = KeyValueField(labels=dict(THREAD_TYPES))
#     privacy = KeyValueField(labels=dict(THREAD_PRIVACIES))
#     link = serializers.HyperlinkedIdentityField(view_name='api:thread-detail')
#
#     class Meta:
#         model = Thread
#         fields = ('id',
#                   'link',
#                   'type',
#                   'title',
#                   'user',
#                   'privacy',
#                   'state',
#                   'paper',
#                   'title',
#                   'content',
#                   'members',
#                   'posts',
#                   'created',
#                   'modified')
#
#     def get_state(self, obj):
#         if ThreadUser.objects.filter(user=self.context['request'].user,
#                                      thread=obj).exists():
#             instance = ThreadUser.objects.get(user=self.context['request'].user,
#                                               thread=obj)
#             return ThreadUserSerializer(instance=instance,
#                                         context={'request': self.context[
#                                             'request']}).data
#         else:
#             return None


class ThreadSerializer(serializers.HyperlinkedModelSerializer):

    state = serializers.SerializerMethodField()
    members = serializers.SerializerMethodField()
    posts = serializers.SerializerMethodField()

    class Meta:
        model = Thread
        extra_kwargs = {
            'link': {'view_name': 'api:thread-detail'},
            'user': {'view_name': 'api:user-detail'},
            'paper': {'view_name': 'api:paper-detail'},
        }
        fields = (
            'link',
            'user',
            'type',
            'privacy',
            'title',
            'content',
            'paper',
            'state',
            'members',
            'posts',
            'created',
            'modified'
        )
        read_only_fields = (
            'user',
            'state',
            'members',
            'posts',
            'created',
            'modified'
        )

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        threaduser = obj.state(self.context['request'].user)
        if threaduser:
            return reverse('api:threaduser-detail',
                           kwargs={'pk': threaduser.id},
                           request=self.context['request'],
                           format=self.context['format'])
        return None

    def get_members(self, obj):
        """Get thread members based on if request.user is himself a Thread member"""
        members = obj.members
        members_urls = []
        if self.context['request'].user in members:
            for member in members:
                members_urls.append(
                    reverse('api:user-detail',
                           kwargs={'pk': member.id},
                           request=self.context['request'],
                           format=self.context['format'])
                )
        else:
            for member in members[:settings.MAX_MEMBERS_NOT_JOINED]:
                members_urls.append(
                    reverse('api:user-detail',
                           kwargs={'pk': member.id},
                           request=self.context['request'],
                           format=self.context['format'])
                )
        return members_urls

    def get_posts(self, obj):
        """Get thread posts based on if request.user is a Thread member"""
        posts_urls = []
        if self.context['request'].user in obj.members:
            for post in obj.posts.all():
                posts_urls.append(
                    reverse('api:threadpost-detail',
                           kwargs={'pk': post.id},
                           request=self.context['request'],
                           format=self.context['format'])
                )
        return posts_urls

    def get_fields(self, *args, **kwargs):
        """Limit paper queryset to user library paper"""
        fields = super(ThreadSerializer, self).get_fields(*args, **kwargs)
        if self.context.get('request'):
            fields['paper'].queryset = self.context['request'].user.lib.papers.all()
        return fields

    def validate(self, data):
        """
        """
        # Check that paper is in user library if <type>='Paper'
        paper = data['paper']
        if dict(THREAD_TYPES).get(data['type']).lower() == 'paper':
            if paper not in self.context['request'].user.lib.papers.all():
                raise serializers.ValidationError("Related paper must be in user library")
        elif dict(THREAD_TYPES).get(data['type']).lower() == 'question':
            if paper is not None:
                raise serializers.ValidationError("Thread type <Question> cannot have related paper")
        return data

    def save(self, **kwargs):
        return super(ThreadSerializer, self).save(user=self.context['request'].user, **kwargs)

