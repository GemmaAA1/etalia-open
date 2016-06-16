# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json

from django.http import JsonResponse
from django.core.urlresolvers import reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings
from django.template.response import TemplateResponse

from .models import Stream
from .tasks import update_stream, update_trend, reset_stream, reset_trend


@login_required
def my_feeds(request):
    return TemplateResponse(
        request,
        'feeds/list.html',
        {'control_states': json.dumps(
            request.session.get('feeds-control-states',
                                {'time-span': settings.FEEDS_DEFAULT_TIME_SPAN,
                                 'search': None,
                                 'pin': 0}))}
    )


@login_required
def reset_stream_view(request, stream_name):
    if request.is_ajax() or settings.DEBUG:
        reset_stream.delay(request.user.pk, stream_name=stream_name,
                           restrict_journal=False)
        data = {'display_update_modal': True,
                'message': 'Stream reset launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_stream_view(request, stream_name):
    if request.is_ajax() or settings.DEBUG:
        update_stream.delay(request.user.pk, stream_name=stream_name)
        data = {'display_update_modal': True,
                'message': 'Stream update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_trend_view(request, trend_name):
    if request.is_ajax() or settings.DEBUG:
        update_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def reset_trend_view(request, trend_name):
    if request.is_ajax() or settings.DEBUG:
        reset_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def ajax_user_feed_message(request, name):
    if request.method == 'GET':
        userfeed = get_object_or_404(Stream,
                                     name=name,
                                     user=request.user)
        if userfeed.state == 'IDL':
            data = {'done': True,
                    'url': str(reverse_lazy('feeds:feed',
                                        kwargs={'name': userfeed.name}))}
        else:
            data = {'done': False,
                    'message': userfeed.message}
        return JsonResponse(data)



