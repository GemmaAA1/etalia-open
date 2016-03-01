# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import template

register = template.Library()


@register.inclusion_tag('_fragment/loading.html')
def show_loading(type):
    return {'type': type}
