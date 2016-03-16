# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
from django.conf import settings
from django.core.urlresolvers import reverse


def admin_context(request):
    # return the value you want as a dictionnary. you may add multiple values in there.
    return {'environment': os.path.splitext(os.path.basename(settings.CONFIG_FILE))[0],
            'hide_cluster_icon': settings.HIDE_CLUSTER_ICON,
            'version': settings.VERSION
            }


def user_update_check(request):
    """Check for stream, trend update and user initialization"""

    STREAM_CHECK_URL = [
        reverse('feeds:stream'),
        reverse('user:settings'),
    ]

    TREND_CHECK_URL = [
        reverse('feeds:trend'),
        reverse('user:settings'),
    ]

    if not request.user.is_anonymous():
        if not request.user.init_step == 'IDL':
            return {'busy_url': reverse('user:user-init-check')}
        elif request.path in STREAM_CHECK_URL and \
                not request.user.streams.first().state == 'IDL':
            return {'busy_url': reverse('user:user-update-stream-check')}
        elif request.path in TREND_CHECK_URL and \
                not request.user.trends.first().state == 'IDL':
            return {'busy_url': reverse('user:user-update-trend-check')}

    return {'busy_url': None}
