# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.template.response import TemplateResponse


def my_papers(request):
    return TemplateResponse(
        request,
        'papers/my_list.html',
        {'control_states': json.dumps(
            request.session.get('library-control-states',
                                {'time-span': None,
                                 'search': None,
                                 'pin': 0}
                                )
        )
        }
    )


def papers(request):
    return TemplateResponse(
        request,
        'papers/list.html',
        {}
    )