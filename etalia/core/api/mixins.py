# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class MultiSerializerMixin(object):

    exclude_action_serializers = {}

    def get_serializer_class(self):
        view = self.request.query_params.get('view', None)

        if view:
            exclude_serializer = \
                self.exclude_action_serializers.get(self.action, [])
            if view not in exclude_serializer:
                return self.serializer_class.get(view,
                                                 self.serializer_class['default'])

        return self.serializer_class.get(self.action,
                                         self.serializer_class['default'])



