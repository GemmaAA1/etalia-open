# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

READ_METHODS = ('GET', 'HEAD', 'OPTIONS')
WRITE_METHODS = ('PATCH', 'POST', 'PUT')


class MultiSerializerMixin(object):

    exclude_action_serializers = {}
    serializer_class = {}
    request = None
    action = ''

    def get_serializer_class(self):
        # Based of view arguments
        view = self.request.query_params.get('view', None)
        if view:
            exclude_serializer = \
                self.exclude_action_serializers.get(self.action, [])
            if view not in exclude_serializer:
                return self.serializer_class.get(view,
                                                 self.serializer_class['default'])

        if self.serializer_class.get(self.action, None):
            return self.serializer_class.get(self.action)
        return self.serializer_class.get(self.action,
                                         self.serializer_class['default'])


class One2OneNestedLinkSwitchMixin(object):
    """Add support for:
       - [GET, OPTION, HEAD]: one2one relationships are nested by default
       - [POST, PATCH, PUT]: support both, HyperlinkedRelated URL or nested
       representation
       """

    one2one_nested = True

    def __init__(self, *args, **kwargs):
        if 'one2one_nested' in kwargs:
            self.one2one_nested = kwargs.pop('one2one_nested')
        meta = getattr(self, 'Meta', None)
        if not getattr(meta, 'switch_kwargs', None):
            meta.switch_kwargs = {}
        super(One2OneNestedLinkSwitchMixin, self).__init__(*args, **kwargs)

    def get_fields(self, *args, **kwargs):
        """Limit paper queryset to user library paper"""
        fields = super(One2OneNestedLinkSwitchMixin, self).get_fields(*args, **kwargs)
        method = self.context['request'].method
        if self.one2one_nested:
            for field in self.Meta.switch_kwargs.keys():
                clean_kwargs = self.Meta.switch_kwargs[field].copy()
                serializer = clean_kwargs.pop('serializer', None)
                clean_kwargs.pop('queryset', None)
                if method in READ_METHODS:
                    fields[field] = serializer(many=False, read_only=True, **clean_kwargs)
                elif method in WRITE_METHODS and hasattr(self, 'initial_data'):
                    if field in self.initial_data:
                        value = self.initial_data.get(field, None)
                        # nested
                        if isinstance(value, dict):
                            fields[field] = serializer(many=False, read_only=True,
                                                       **clean_kwargs)
                # add queryset
                if 'queryset' in self.Meta.switch_kwargs[field]:
                    queryset = self.get_field_queryset(field)
                    fields[field].queryset = queryset

        return fields

    def get_field_queryset(self, field):
        if 'queryset' in self.Meta.switch_kwargs[field]:
            if isinstance(self.Meta.switch_kwargs[field].get('queryset'), str):
                # method based
                queryset = eval('self.{func_name}()'.format(
                    func_name=self.Meta.switch_kwargs[field].get('queryset')))
                return queryset
            else:
                return self.Meta.switch_kwargs[field].get('queryset')
        return None




