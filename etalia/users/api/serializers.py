# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from django.contrib.auth import get_user_model
from ..models import UserLibPaper, UserLib, Relationship
from etalia.library.api.serializers import PaperNestedSerializer

User = get_user_model()


class UserSerializer(serializers.HyperlinkedModelSerializer):
    photo_url = serializers.URLField(read_only=True, source='photo.name')

    class Meta:
        model = User
        extra_kwargs = {
            'link': {'view_name': 'api:user-detail'},
        }
        fields = (
            'id',
            'link',
            'email',
            'first_name',
            'last_name',
            'photo_url')
        read_only_fields = (
            '__all__'
        )


class UserFullSerializer(serializers.HyperlinkedModelSerializer):
    photo_url = serializers.URLField(read_only=True, source='photo.name')
    user_lib = serializers.HyperlinkedRelatedField(
        read_only=True,
        view_name='api:userlib-detail',
        source='lib')
    following = serializers.SerializerMethodField()
    followers = serializers.SerializerMethodField()

    class Meta:
        model = User
        extra_kwargs = {
            'link': {'view_name': 'api:user-detail'},
            'user_lib': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'id',
            'link',
            'email',
            'first_name',
            'last_name',
            'photo_url',
            'user_lib',
            'following',
            'followers')
        read_only_fields = (
            '__all__'
        )

    def get_following(self, obj):
        following = obj.following
        following_list = []
        for user in following:
            following_list.append(
                UserSerializer(
                    instance=user,
                    context={'request': self.context['request']}).data
            )
        return following_list

    def get_followers(self, obj):
        followers = obj.followers
        followers_list = []
        for user in followers:
            followers_list.append(
                UserSerializer(
                    instance=user,
                    context={'request': self.context['request']}).data
            )
        return followers_list


class UserLibPaperSerializer(serializers.HyperlinkedModelSerializer):

    user_lib = serializers.HyperlinkedRelatedField(source='userlib',
                                                   view_name='api:userlib-detail',
                                                   read_only=True)

    class Meta:
        model = UserLibPaper
        extra_kwargs = {
            'link': {'view_name': 'api:userlibpaper-detail'},
            'paper': {'view_name': 'api:paper-detail'},
        }
        fields = (
            'id',
            'link',
            'date_created',
            'authored',
            'user_lib',
            'paper')
        read_only_fields = (
            '__all__'
        )


class UserLibPaperNestedSerializer(serializers.HyperlinkedModelSerializer):

    paper = PaperNestedSerializer(many=False, read_only=True)

    user_lib = serializers.HyperlinkedRelatedField(source='userlib',
                                                   view_name='api:userlib-detail',
                                                   read_only=True)

    class Meta:
        model = UserLibPaper
        extra_kwargs = {
            'link': {'view_name': 'api:userlibpaper-detail'},
        }
        fields = (
            'id',
            'link',
            'date_created',
            'authored',
            'user_lib',
            'paper')
        read_only_fields = (
            '__all__'
        )


class UserLibSerializer(serializers.HyperlinkedModelSerializer):

    user = serializers.HyperlinkedRelatedField(view_name='api:user-detail',
                                               read_only=True)
    user_lib_papers = UserLibPaperSerializer(many=True, read_only=True,
                                             source='userlib_paper')
    id = serializers.IntegerField(source='user_id')

    class Meta:
        model = UserLib
        extra_kwargs = {
            'link': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'user_lib_papers',
        )
        read_only_fields = (
            '__all__'
        )


class UserLibNestedSerializer(serializers.HyperlinkedModelSerializer):

    user = UserSerializer(many=False, read_only=True)
    user_lib_papers = UserLibPaperNestedSerializer(many=True, read_only=True,
                                                   source='userlib_paper')
    id = serializers.IntegerField(source='user_id')

    class Meta:
        model = UserLib
        extra_kwargs = {
            'link': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'user_lib_paper',
        )
        read_only_fields = (
            '__all__'
        )


class RelationshipSerializer(serializers.HyperlinkedModelSerializer):

    # status

    class Meta:
        model = Relationship
        extra_kwargs = {
            'link': {'view_name': 'api:relationship-detail'},
            'from_user': {'view_name': 'api:user-detail'},
            'to_user': {'view_name': 'api:user-detail'}
        }
        fields = (
            'id',
            'link',
            'from_user',
            'to_user',
            'status'
        )

    def validate_from_user(self, value):
        """From_user must be the request.user"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "from_user deserialized and user logged in are different")
        return value

    def validate_to_user(self, value):
        """ Cannot have a relationship with yourself, that would be awkward
        """
        if value == self.context['request'].user:
            raise serializers.ValidationError("to_user cannot be request.user")
        return value
