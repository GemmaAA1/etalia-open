# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
from django.conf import settings
from django.core.urlresolvers import reverse


def admin_context(request):
    # return the value you want as a dictionnary. you may add multiple values in there.
    return {'environment': os.path.splitext(os.path.basename(settings.CONFIG_FILE))[0],
            'hide_cluster_icon': settings.HIDE_CLUSTER_ICON}


def user_update_check(request):
    """Check for stream, trend update and user initialization"""

    if not request.user.init_step == 'IDL':
        return {'busy_url': reverse('user:user-init-check')}
    elif request.path.startswith('/feed/stream') and \
            not request.user.streams.first().state == 'IDL':
        return {'busy_url': reverse('user:user-update-stream-check')}
    elif request.path.startswith('/feed/trend') and \
            not request.user.trends.first().state == 'IDL':
        return {'busy_url': reverse('user:user-update-trend-check')}
    else:
        return {'busy_url': None}
