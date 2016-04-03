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
    lib = serializers.HyperlinkedRelatedField(read_only=True,
                                              view_name='api:userlib-detail')

    class Meta:
        model = User
        extra_kwargs = {
            'link': {'view_name': 'api:user-detail'},
            'lib': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'id',
            'link',
            'email',
            'first_name',
            'last_name',
            'photo_url',
            'lib')
        read_only_fields = (
            '__all__'
        )



class UserLibPaperSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = UserLibPaper
        extra_kwargs = {
            'link': {'view_name': 'api:userlibpaper-detail'},
            'paper': {'view_name': 'api:paper-detail'},
            'userlib': {'view_name': 'api:userlib-detail'}
        }
        fields = (
            'id',
            'link',
            'date_created',
            'authored',
            'userlib',
            'paper')
        read_only_fields = (
            '__all__'
        )


class UserLibPaperNestedSerializer(serializers.HyperlinkedModelSerializer):

    paper = PaperNestedSerializer(many=False, read_only=True)

    class Meta:
        model = UserLibPaper
        extra_kwargs = {
            'link': {'view_name': 'api:userlibpaper-detail'},
            'userlib': {'view_name': 'api:userlib-detail'}
        }
        fields = (
            'id',
            'link',
            'date_created',
            'authored',
            'userlib',
            'paper')
        read_only_fields = (
            '__all__'
        )


class UserLibSerializer(serializers.HyperlinkedModelSerializer):

    user = serializers.HyperlinkedRelatedField(view_name='api:user-detail',
                                               read_only=True)
    userlib_paper = UserLibPaperSerializer(many=True, read_only=True)

    class Meta:
        model = UserLib
        extra_kwargs = {
            'link': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'link',
            'user',
            'userlib_paper',
        )
        read_only_fields = (
            '__all__'
        )


class UserLibNestedSerializer(serializers.HyperlinkedModelSerializer):

    user = UserSerializer(many=False, read_only=True)
    userlib_paper = UserLibPaperNestedSerializer(many=True, read_only=True)

    class Meta:
        model = UserLib
        extra_kwargs = {
            'link': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'link',
            'user',
            'userlib_paper',
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
