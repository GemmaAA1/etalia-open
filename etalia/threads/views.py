# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.template.response import TemplateResponse


def thread(request, pk):
    return TemplateResponse(
        request,
        'threads/detail.html',
        {'id': pk}
    )


def my_threads(request):
    return TemplateResponse(
        request,
        'threads/my_list.html',
        {'control_states': json.dumps(request.session.get(
            'threads-control-states',
            {'time-span': None, 'search': None, 'pin': 0}))}
    )


def threads(request):
    return TemplateResponse(
        request,
        'threads/list.html',
        {}
    )
