# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.template.defaulttags import register

@register.filter
def get_is_liked(dictionary, key):
    if dictionary.get(key, None):
        return dictionary.get(key).get('liked')
    else:
        return False

@register.filter
def get_is_disliked(dictionary, key):
    if dictionary.get(key, None):
        return dictionary.get(key).get('disliked')
    else:
        return False
