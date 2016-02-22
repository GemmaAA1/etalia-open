# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
from django.conf import settings


def admin_context(request):
    # return the value you want as a dictionnary. you may add multiple values in there.
    return {'environment': os.path.splitext(os.path.basename(settings.CONFIG_FILE))[0]}
