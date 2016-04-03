# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf import settings

from rest_framework import serializers
from rest_framework.reverse import reverse
from rest_framework.fields import empty

from etalia.users.api.serializers import UserSerializer
from etalia.library.api.serializers import PaperSerializer, PaperNestedSerializer
from ..models import Thread, ThreadPost, ThreadComment, ThreadUser
from ..constant import THREAD_PRIVACIES, THREAD_TYPES


class PatchSerializer(serializers.Serializer):

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


class ThreadUserSerializer(serializers.HyperlinkedModelSerializer):

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

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged are different")
        return value


class ThreadCommentSerializer(serializers.HyperlinkedModelSerializer):
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

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user loged in are different")
        return value

    def validate(self, data):
        """User is a member of thread"""
        thread = data['thread']
        if not self.context['request'].user in thread.members:
            raise serializers.ValidationError(
                "user is not a member of this thread")
        return data


class ThreadCommentNestedSerializer(serializers.HyperlinkedModelSerializer):
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


class ThreadPostSerializer(serializers.HyperlinkedModelSerializer):
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
            'modified'
        )
        read_only_fields = (
            'id',
            'link',
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

    def get_fields(self, *args, **kwargs):
        """Limit paper queryset to user library paper"""
        fields = super(ThreadSerializer, self).get_fields(*args, **kwargs)
        if self.context.get('request'):
            fields['paper'].queryset = self.context[
                'request'].user.lib.papers.all()
        return fields

    def validate_user(self, value):
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged in are different")
        return value

    def validate(self, data):
        """Integrity of thread <type> and <paper>
        """
        # Check that paper is in user library if <type>='Paper'
        paper = data['paper']
        if dict(THREAD_TYPES).get(data['type']).lower() == 'paper':
            if paper not in self.context['request'].user.lib.papers.all():
                raise serializers.ValidationError(
                    "Related paper must be in user library")
        elif dict(THREAD_TYPES).get(data['type']).lower() == 'question':
            if paper is not None:
                raise serializers.ValidationError(
                    "Thread type <Question> cannot have related paper")
        return data


class ThreadNestedSerializer(serializers.HyperlinkedModelSerializer):
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
            'modified'
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
