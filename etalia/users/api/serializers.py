# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from django.contrib.auth import get_user_model
from ..models import UserLibPaper, UserLib, Relationship
from etalia.library.api.serializers import PaperSerializer

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

class UserLibPaperSerializer(serializers.HyperlinkedModelSerializer):
    paper = PaperSerializer(many=False, read_only=True)

    class Meta:
        model = UserLibPaper
        extra_kwargs = {
            'link': {'view_name': 'api:userlib-detail'},
        }
        fields = (
            'id',
            'date_created',
            'authored',
            'paper')
        read_only_fields = (
            '__all__'
        )


class UserLibSerializer(serializers.ModelSerializer):
    papers = UserLibPaperSerializer(many=True, read_only=True,
                                    source='userlib_paper')

    class Meta:
        model = UserLib
        fields = ('pk',
                  'papers')
        read_only_fields = (
            '__all__'
        )


class RelationshipSerializer(serializers.ModelSerializer):
    class Meta:
        model = Relationship
        fields = ('pk',
                  'from_user',
                  'to_user')
