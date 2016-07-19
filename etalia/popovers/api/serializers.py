# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from django.template.loader import render_to_string, TemplateDoesNotExist
from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin

from ..models import PopOver, UserPopOver, UserPopOverUpdateDisplay
from ..constants import POPOVER_TYPES, ANCHORED, MODAL, GOT_IT


class PopOverSerializer(serializers.HyperlinkedModelSerializer):

    body = serializers.SerializerMethodField()

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
            'priority',
            'template_path',
        )
        read_only_fields = (
            'body',
        )

    def get_body(self, obj):
        """Render template file"""
        # template = get_template(obj.template_path)
        context = {}
        try:
            return render_to_string(obj.template_path, context,
                                request=self.context['request'])
        except TemplateDoesNotExist:
            return 'No template file match'


    def validate(self, data):
        """Check that anchored is defined is type is Anchored and is empty if type is Modal"""
        if 'type' in data:
            type = data.get('type', None)
            anchor = data.get('anchor', None)
            if type == ANCHORED:
                if not anchor:
                    raise serializers.ValidationError(
                        "Anchor must be defined with PopOver of type Anchored")
            elif type == MODAL:
                if anchor:
                    data['anchor'] = ''
        return data


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
            'status',
            'display',
        )
        read_only_fields = (
            'id',
            'link',
            'user',
            'popover',
            'display',
        )
        switch_kwargs = {
            'popover': {'serializer': PopOverSerializer},
        }

    def save(self, **kwargs):
        super(UserPopOverSerializer, self).save(**kwargs)
        # If user GOT_IT, triggered deferred update of display popovers
        if self.validated_data.get('status') == GOT_IT:
            upoud, _ = UserPopOverUpdateDisplay.objects\
                .get_or_create(user=self.instance.user)
            upoud.deferred_display_update()

        return self.instance
