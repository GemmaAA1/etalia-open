# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf import settings
from django.contrib.auth import get_user_model

from rest_framework import serializers
from rest_framework.reverse import reverse
from rest_framework.fields import empty

from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin
from etalia.users.api.serializers import UserSerializer, UserFilterSerializer
from etalia.library.api.serializers import PaperSerializer, PaperNestedSerializer
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser, ThreadUserInvite
from ..constant import THREAD_PRIVACIES, THREAD_TYPES, THREAD_PRIVATE

User = get_user_model()

WRITE_METHODS = ('PATCH', 'POST', 'PUT')
READ_METHODS = ('GET', 'HEAD', 'OPTIONS')


class PatchSerializer(serializers.Serializer):
    """Serializer for PATCH update based on operations"""
    OPERATIONS = (
        "replace",
    )

    op = serializers.ChoiceField(choices=OPERATIONS, required=True)
    path = serializers.URLField(required=False)
    value = serializers.IntegerField(required=False)

    class Meta:
        fields = (
            'op',
            'path',
            'value',
        )


class ThreadFilterSerializer(serializers.BaseSerializer):
    """Serializer for filters on side panel of Threads list"""
    def to_representation(self, instance):
        return {
            'users': [UserFilterSerializer(instance=user, context=self.context).data
                      for user in instance.get('users', None)]
        }


class ThreadUserInviteSerializer(One2OneNestedLinkSwitchMixin,
                                 serializers.HyperlinkedRelatedField):
    """ThreadUserInvite serializer"""

    class Meta:
        model = ThreadUserInvite
        extra_kwargs = {

        }



class ThreadUserSerializer(One2OneNestedLinkSwitchMixin,
                           serializers.HyperlinkedModelSerializer):
    """ThreadUser serializer"""
    class Meta:
        model = ThreadUser
        extra_kwargs = {
            'link': {'view_name': 'api:threaduser-detail'},
            'user': {'view_name': 'api:user-detail'},
            'thread': {'view_name': 'api:thread-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'thread',
            'watch',
            'participate'
        )
        read_only_fields = (
            'id',
            'link',
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer},
        }

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged are different")
        return value

    def validate_thread(self, value):
        if not value.published_at:
            raise serializers.ValidationError("Thread is not yet published")
        if value.privacy == THREAD_PRIVATE:
            # User must have received an invite to interact with thread or be the owner
            if not ThreadUserInvite.objects.exists(thread=value,
                                                   to_user=self.context['request'].user) \
                    or value.user == self.context['request'].user:
                raise serializers.ValidationError("Thread is private. You need an invite to join")
        return value


class ThreadCommentSerializer(One2OneNestedLinkSwitchMixin,
                              serializers.HyperlinkedModelSerializer):
    """ThreadComment serializer"""
    class Meta:
        model = ThreadComment
        extra_kwargs = {
            'link': {'view_name': 'api:threadcomment-detail'},
            'user': {'view_name': 'api:user-detail'},
            'post': {'view_name': 'api:threadpost-detail'}
        }
        fields = (
            'id',
            'link',
            'user',
            'post',
            'content',
            'created',
            'modified',
        )
        read_only_fields = (
            'id',
            'link',
            'created',
            'modified')

        switch_kwargs = {
            'user': {'serializer': UserSerializer},
        }

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user loged in are different")
        return value

    def validate(self, data):
        """User is a member of thread"""
        post = data['post']
        if not self.context['request'].user in post.thread.members:
            raise serializers.ValidationError(
                "user is not a member of this thread")
        return data


class ThreadCommentNestedSerializer(serializers.HyperlinkedModelSerializer):
    """ThreadComment nested serializer"""

    user = UserSerializer(many=False, read_only=True)

    class Meta:
        model = ThreadComment
        extra_kwargs = {
            'link': {'view_name': 'api:threadcomment-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'content',
            'created',
            'modified',
        )
        read_only_fields = (
            '__all__'
        )


class ThreadPostSerializer(One2OneNestedLinkSwitchMixin,
                           serializers.HyperlinkedModelSerializer):
    """ThreadPost serializer"""
    class Meta:
        model = ThreadPost
        extra_kwargs = {
            'link': {'view_name': 'api:threadpost-detail'},
            'comments': {'view_name': 'api:threadcomment-detail'},
            'user': {'view_name': 'api:user-detail', 'required': True},
            'thread': {'view_name': 'api:thread-detail'}
        }
        fields = (
            'id',
            'link',
            'thread',
            'user',
            'content',
            'created',
            'modified',
            'comments')
        read_only_fields = (
            'id',
            'link',
            'created',
            'modified',
            'comments')

        switch_kwargs = {
            'user': {'serializer': UserSerializer}
        }

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged in are different")
        return value

    def validate(self, data):
        """User is a member of thread"""
        thread = data['thread']
        if not self.context['request'].user in thread.members:
            raise serializers.ValidationError(
                "user is not a member of this thread")
        return data


class ThreadPostNestedSerializer(serializers.HyperlinkedModelSerializer):
    """ThreadPost nested serializer"""

    comments = ThreadCommentNestedSerializer(many=True, read_only=True)
    user = UserSerializer(many=False, read_only=True)

    class Meta:
        model = ThreadPost
        extra_kwargs = {
            'link': {'view_name': 'api:threadpost-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'content',
            'created',
            'modified',
            'comments')
        read_only_fields = (
            '__all__')


class ThreadSerializer(One2OneNestedLinkSwitchMixin,
                       serializers.HyperlinkedModelSerializer):
    """Thread serializer"""

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
            'id',
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
            'modified',
            'published_at'
        )
        read_only_fields = (
            'id',
            'link',
            'state',
            'members',
            'posts',
            'created',
            'modified',
            'published_at'
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer},
            'paper': {'serializer': PaperNestedSerializer,
                      'queryset': 'get_paper_queryset',
                      'allow_null': True}
        }

    def get_paper_queryset(self):
        return self.context['request'].user.lib.papers.all()

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        threaduser = obj.state(self.context['request'].user)
        if threaduser:
            return ThreadUserSerializer(
                instance=threaduser,
                context={'request': self.context['request']}).data
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
        """Retrieve posts if request.user is in Thread member"""
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

    def validate_user(self, value):
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged in are different")
        return value

    def validate(self, data):
        """Integrity of thread <type> and <paper>
        """
        # Check that paper is in user library if <type>='Paper'
        if 'paper' in data:
            paper = data.get('paper', None)
            if dict(THREAD_TYPES).get(data.get('type', '')).lower() == 'paper':
                if paper not in self.context['request'].user.lib.papers.all():
                    raise serializers.ValidationError(
                        "Related paper must be in user library")
            elif dict(THREAD_TYPES).get(data.get('type', '')).lower() == 'question':
                if paper is not None:
                    raise serializers.ValidationError(
                        "Thread type <Question> cannot have related paper")
        return data


class ThreadNestedSerializer(serializers.HyperlinkedModelSerializer):
    """Thread nested serializer"""

    state = serializers.SerializerMethodField()
    members = serializers.SerializerMethodField()
    posts = serializers.SerializerMethodField()
    user = UserSerializer(many=False, read_only=True)
    paper = PaperNestedSerializer(many=False, read_only=True)

    class Meta:
        model = Thread
        extra_kwargs = {
            'link': {'view_name': 'api:thread-detail'},
            'user': {'view_name': 'api:user-detail'},
            'paper': {'view_name': 'api:paper-detail'},
            'posts': {'view_name': 'api:threadpost-detail'},
        }
        fields = (
            'id',
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
            'modified',
            'published_at'
        )
        read_only_fields = (
            '__all__',
        )

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        threaduser = obj.state(self.context['request'].user)
        if threaduser:
            return ThreadUserSerializer(
                instance=threaduser,
                context={'request': self.context['request']}).data
        return None

    def get_members(self, obj):
        """Get thread members based on if request.user is himself a Thread member"""
        members = obj.members
        members_list = []
        if self.context['request'].user in members:
            for member in members:
                members_list.append(
                    UserSerializer(
                        instance=member,
                        context={'request': self.context['request']}).data
                )
        else:
            for member in members[:settings.MAX_MEMBERS_NOT_JOINED]:
                members_list.append(
                    UserSerializer(
                        instance=member,
                        context={'request': self.context['request']}).data
                )
        return members_list

    def get_posts(self, obj):
        """Retrieve posts if request.user is in Thread member"""
        posts_list = []
        if self.context['request'].user in obj.members:
            for post in obj.posts.all():
                posts_list.append(
                    ThreadPostNestedSerializer(
                        instance=post,
                        context={'request': self.context['request']}).data
                )
        return posts_list
