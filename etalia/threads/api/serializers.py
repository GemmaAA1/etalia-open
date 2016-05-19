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
from ..constant import THREAD_PRIVACIES, THREAD_TYPES, THREAD_PRIVATE, \
    THREAD_JOINED, THREAD_LEFT, THREAD_INVITE_PENDING, THREAD_INVITE_ACCEPTED

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
            'groups': [
                {
                    "name": "user_id",
                    "label": "Users",
                    "entries": [
                        UserFilterSerializer(instance=user,
                                             context=self.context).data
                        for user in instance.get('users', None)
                        ]
                }
            ]
        }


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
                      'allow_null': True,
                      }
        }

    def get_paper_queryset(self):
        return self.context['request'].user.lib.papers.all()

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        threaduser = obj.state(self.context['request'].user)
        if threaduser:
            if self.one2one_nested:
                return ThreadUserSerializer(
                    instance=threaduser,
                    context={'request': self.context['request']},
                    one2one_nested=False
                ).data
            else:
                return reverse('api:threaduser-detail',
                               kwargs={'pk': threaduser.id},
                               request=self.context.get('request', None),
                               format=self.context.get('format', None))
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
                            request=self.context.get('request', None),
                            format=self.context.get('format', None))
                )
        else:
            for member in members[:settings.MAX_MEMBERS_NOT_JOINED]:
                members_urls.append(
                    reverse('api:user-detail',
                            kwargs={'pk': member.id},
                            request=self.context.get('request', None),
                            format=self.context.get('format', None))
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
                            request=self.context.get('request', None),
                            format=self.context.get('format', None))
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
            'user': {'serializer': UserSerializer,
                     'one2one_nested': False},
            'thread': {'serializer': ThreadSerializer,
                       'one2one_nested': False},
        }

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged are different")
        return value

    def validate_participate(self, value):
        """Check that thread is joined if request is to leave thread"""
        if self.initial_data:
            if value == THREAD_LEFT and \
                    not self.initial_data.get('participate') == THREAD_JOINED:
                raise serializers.ValidationError(
                    "You cannot leave a thread that you have not joined")
        return value

    def validate_thread(self, value):
        if not value.published_at:
            raise serializers.ValidationError("Thread is not yet published")
        if value.privacy == THREAD_PRIVATE:
            # User must have received an invite to interact with thread or be the owner
            if not ThreadUserInvite.objects.filter(thread=value,
                                                   to_user=self.context['request'].user).exists() \
                    or value.user == self.context['request'].user:
                raise serializers.ValidationError("Thread is private. You need an invite to join")
        return value


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
            'user': {'serializer': UserSerializer},
            'thread': {'serializer': ThreadSerializer,
                       'one2one_nested': False}
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
            'post': {'serializer': ThreadPostSerializer,
                     'one2one_nested': False}
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


class ThreadPostNestedSerializer(ThreadPostSerializer):
    """ThreadPost nested serializer"""

    comments = ThreadCommentSerializer(many=True, read_only=True)
    user = UserSerializer(many=False, read_only=True)

    class Meta(ThreadPostSerializer.Meta):
        read_only_fields = (
            '__all__')


class ThreadNestedSerializer(ThreadSerializer):
    """Thread nested serializer"""

    state = serializers.SerializerMethodField()
    members = serializers.SerializerMethodField()
    posts = serializers.SerializerMethodField()
    paper = PaperNestedSerializer(many=False, read_only=True)

    class Meta(ThreadSerializer.Meta):
        read_only_fields = (
            '__all__',
        )

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


class ThreadUserInviteSerializer(One2OneNestedLinkSwitchMixin,
                                 serializers.HyperlinkedModelSerializer):
    """ThreadUserInvite serializer"""

    class Meta:
        model = ThreadUserInvite
        extra_kwargs = {
            'link': {'view_name': 'api:threaduserinvite-detail'},
            'to_user': {'view_name': 'api:user-detail'},
            'from_user': {'view_name': 'api:user-detail'},
            'thread': {'view_name': 'api:thread-detail'},
        }
        fields = (
            'link',
            'id',
            'thread',
            'from_user',
            'to_user',
            'status',
        )
        switch_kwargs = {
            'to_user': {'serializer': UserSerializer},
            'from_user': {'serializer': UserSerializer},
            'thread': {'serializer': ThreadSerializer,
                       'one2one_nested': False},
        }
        read_only_fields = (
            'id',
            'link',
        )
        validators = []

    def validate_thread(self, value):
        # If thread private, must be owner to invite people to thread
        if value.privacy == THREAD_PRIVATE and \
                        not value.user == self.context['request'].user:
            raise serializers.ValidationError(
                "You cannot invite another user to join this thread. "
                "This thread is private and your are not the owner")
        # Thread must be published
        if value.published_at is None:
            raise serializers.ValidationError(
                "You cannot invite another user to this thread. "
                "The thread is not yet published.")
        # User must be a member of thread
        if self.context['request'].user not in value.members:
            raise serializers.ValidationError(
                "You cannot invite another user to this thread. "
                "Your are not a member of this thread")

        return value

    def validate_from_user(self, value):
        """from_user must be the current user"""
        if not self.initial_data and not value == self.context['request'].user:
            raise serializers.ValidationError(
                "from_user deserialized and current user are different")
        return value

    def validate_to_user(self, value):
        """ Cannot have an invite to yourself or to a member of the thread"""
        if self.context['view'].action == 'create':
            if value == self.context['request'].user:
                raise serializers.ValidationError("to_user cannot be current user")
        return value

    def validate_status(self, value):
        """ Status is PENDING for new invite"""
        if self.context['view'].action == 'create':
            if not value == THREAD_INVITE_PENDING:
                raise serializers.ValidationError("status must be PENDING for new invite")
        return value

    def validate(self, data):
        """Validate that:
            user has not already been invited
            user is not a member of the thread already"""
        if self.context['view'].action == 'create':
            # user has not already been invited
            if ThreadUserInvite.objects.filter(thread=data.get('thread'),
                                               from_user=data.get('from_user'),
                                               to_user=data.get('to_user'))\
                    .exists():
                raise serializers.ValidationError("You already invited this user to this thread")
            # user is not a member of the thread
            to_user = data.get('to_user')
            thread = data.get('thread')
            if to_user in thread.members:
                raise serializers.ValidationError("The user you are inviting is already in the thread")
        return data

    def save(self, **kwargs):
        instance = super(ThreadUserInviteSerializer, self).save(**kwargs)
        # When accepted, set related ThreadUser to joined
        if instance.status == THREAD_INVITE_ACCEPTED:
            tu, _ = ThreadUser.objects.get_or_create(user=instance.to_user,
                                                     thread=instance.thread)
            tu.participate = THREAD_JOINED
            tu.save()


class ThreadUserInviteUpdateSerializer(ThreadUserInviteSerializer):
    """ThreadUserInvite serializer"""

    class Meta(ThreadUserInviteSerializer.Meta):
        read_only_fields = (
            'id',
            'link',
            'thread',
            'from_user',
            'to_user'
        )
