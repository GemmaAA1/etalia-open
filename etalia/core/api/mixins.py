# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class MultiSerializerMixin(object):

    def get_serializer_class(self):
        view = self.request.query_params.get('view', None)
        if view:
            return self.serializer_class.get(view,
                                             self.serializer_class['default'])
        else:
            return self.serializer_class.get(self.action,
                                             self.serializer_class['default'])



