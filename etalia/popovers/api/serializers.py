# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin

from ..models import PopOver, UserPopOver


class PopOverSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:

        model = PopOver
        extra_kwargs = {
            'link': {'view_name': 'api:popover-detail'}
        }
        fields = (
            'id',
            'link',
            'title',
            'body',
            'anchor',
            'type',
        )
        # read_only_fields = (
        #     '__all__',
        # )


class UserPopOverSerializer(One2OneNestedLinkSwitchMixin,
                            serializers.HyperlinkedModelSerializer):

    class Meta:

        model = UserPopOver
        extra_kwargs = {
            'link': {'view_name': 'api:userpopover-detail'},
            'popover': {'view_name': 'api:popover-detail'},
            'user': {'view_name': 'api:user-detail'}
        }
        fields = (
            'id',
            'link',
            'user',
            'popover',
            'status'
        )
        read_only_fields = (
            'id',
            'link',
            'user',
            'popover',
        )
        switch_kwargs = {
            'popover': {'serializer': PopOverSerializer},
        }
