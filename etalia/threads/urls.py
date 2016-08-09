# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views


urlpatterns = [
    url(r'^my-threads$', views.my_threads, name='my_threads'),
    url(r'^threads$', views.threads, name='threads'),
    url(r'^(?P<pk>[0-9]+)/$', views.thread, name='thread'),
]