# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.template.defaulttags import register

@register.filter
def get_is_pinned(dictionary, key):
    if dictionary.get(key, None):
        return dictionary.get(key).get('is_pinned')
    else:
        return False

@register.filter
def get_is_banned(dictionary, key):
    if dictionary.get(key, None):
        return dictionary.get(key).get('is_banned')
    else:
        return False

@register.filter
def get_is_in_lib(dictionary, key):
    if dictionary.get(key, None):
        return dictionary.get(key).get('in_lib')
    else:
        return False