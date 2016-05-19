# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from rest_framework.reverse import reverse

from django.contrib.auth import get_user_model

from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin
from etalia.library.api.serializers import PaperNestedSerializer
from ..models import UserLibPaper, UserLib, Relationship
from avatar.models import Avatar

User = get_user_model()


class AvatarSerializer(serializers.ModelSerializer):

    class Meta:
        model = Avatar


class UserSerializer(One2OneNestedLinkSwitchMixin,
                     serializers.HyperlinkedModelSerializer):
    """User serializer"""

    avatars = AvatarSerializer(read_only=True, many=True, source='avatar_set')

    class Meta:
        model = User
        extra_kwargs = {
            'link': {'view_name': 'api:user-detail'},
            'avatar': {'view_name': 'api:avatar-detail'},
        }
        fields = (
            'id',
            'link',
            'email',
            'first_name',
            'last_name',
            'avatar')
        read_only_fields = (
            '__all__'
        )


class UserFilterSerializer(serializers.ModelSerializer):
    """User serializer user in filter side panel"""

    count = serializers.IntegerField(read_only=True)
    label = serializers.SerializerMethodField()

    class Meta:
        model = User
        fields = (
            'id',
            'label',
            'count')
        read_only_fields = (
            '__all__'
        )

    def get_label(self, obj):
        return '{0} {1}'.format(obj.first_name, obj.last_name)


class UserFullSerializer(serializers.HyperlinkedModelSerializer):
    """User serializer with full representation """

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

    user_lib_papers = serializers.SerializerMethodField()
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

    def get_user_lib_papers(self, obj):
        userlib_paper_urls = []
        for userlib_paper in obj.userlib_paper.all():
            userlib_paper_urls.append(
                reverse('api:userlibpaper-detail',
                        kwargs={'pk': userlib_paper.id},
                        request=self.context.get('request', None),
                        format=self.context.get('format', None))
            )
        return userlib_paper_urls


class UserLibNestedSerializer(UserLibSerializer):
    """UserLib nested serializer"""

    class Meta(UserLibSerializer.Meta):
        pass

    def get_user_lib_papers(self, obj):
        userlib_papers = []
        for userlib_paper in obj.userlib_paper.all():
            userlib_papers.append(
                UserLibPaperSerializer(
                    instance=userlib_paper,
                    context={'request': self.context['request']},
                    one2one_nested=True
                ).data
            )
        return userlib_papers


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
