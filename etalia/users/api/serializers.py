# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from django.contrib.auth import get_user_model

from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin
from etalia.library.api.serializers import PaperNestedSerializer
from ..models import UserLibPaper, UserLib, Relationship

User = get_user_model()


class UserSerializer(serializers.HyperlinkedModelSerializer):
    """User serializer"""

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


class UserFilterSerializer(serializers.HyperlinkedModelSerializer):
    """User serializer user in filter side panel"""

    photo_url = serializers.URLField(read_only=True, source='photo.name')
    count = serializers.IntegerField(read_only=True)
    name = serializers.SerializerMethodField()

    class Meta:
        model = User
        extra_kwargs = {
            'link': {'view_name': 'api:user-detail'},
        }
        fields = (
            'id',
            'link',
            'name',
            'photo_url',
            'count')
        read_only_fields = (
            '__all__'
        )

    def get_name(self, obj):
        return '{0} {1}'.format(obj.first_name, obj.last_name)


class UserFullSerializer(serializers.HyperlinkedModelSerializer):
    """User serializer with full representation """
    photo_url = serializers.URLField(read_only=True, source='photo.name')
    user_lib = serializers.HyperlinkedRelatedField(
        read_only=True,
        view_name='api:userlib-detail',
        source='lib')
    following = UserSerializer(many=True, read_only=True)
    followers = UserSerializer(many=True, read_only=True)

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


class UserLibPaperSerializer(One2OneNestedLinkSwitchMixin,
                             serializers.HyperlinkedModelSerializer):
    """UserLibPaper serializer"""

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
        switch_kwargs = {
            'paper': {'serializer': PaperNestedSerializer}
        }


class UserLibSerializer(One2OneNestedLinkSwitchMixin,
                        serializers.HyperlinkedModelSerializer):
    """UserLib serializer"""

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
        switch_kwargs = {
            'user': {'serializer': UserSerializer}
        }


class UserLibNestedSerializer(One2OneNestedLinkSwitchMixin,
                              serializers.HyperlinkedModelSerializer):
    """UserLib nested serializer"""

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
            'user_lib_paper',
        )
        read_only_fields = (
            '__all__'
        )
        switch_kwargs = {
            'user': {'serializer': UserSerializer}
        }


class RelationshipSerializer(One2OneNestedLinkSwitchMixin,
                             serializers.HyperlinkedModelSerializer):
    """Relationship serializer"""

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
        switch_kwargs = {
            'from_user': {'serializer': UserSerializer},
            'to_user': {'serializer': UserSerializer}
        }

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
