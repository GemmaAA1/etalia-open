# -*- coding: utf-8 -*-
from __future__ import unicode_literals

"""
WSGI config for paperstream project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/howto/deployment/wsgi/
"""

import os
from django.core.wsgi import get_wsgi_application
from dj_static import Cling
# from whitenoise.django import DjangoWhiteNoise

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'paperstream.settings.test')
application = get_wsgi_application()
application = Cling(application)
# application = DjangoWhiteNoise(application)
