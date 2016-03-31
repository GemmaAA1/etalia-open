# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers


class KeyValueField(serializers.Field):
    """ A field that takes a field's value as the key and returns
    the associated value for serialization """

    labels = {}
    inverted_labels = {}

    def __init__(self, labels, *args, **kwargs):
        self.labels = labels
        # Check to make sure the labels dict is reversible, otherwise
        # deserialization may produce unpredictable results
        inverted = {}
        for k, v in labels.items():
            if v in inverted:
                raise ValueError(
                    'The field is not deserializable with the given labels.'
                    ' Please ensure that labels map 1:1 with values'
                )
            inverted[v] = k
        self.inverted_labels = inverted
        return super(KeyValueField, self).__init__(*args, **kwargs)

    def to_representation(self, obj):
        if type(obj) is list:
            return [self.labels.get(o, None) for o in obj]
        else:
            return self.labels.get(obj, None)

    def to_internal_value(self, data):
        if type(data) is list:
            return [self.inverted_labels.get(o, None) for o in data]
        else:
            return self.inverted_labels.get(data, None)
