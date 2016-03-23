# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views


urlpatterns = [
    url(r'^create$', views.create, name='create'),
    url(r'^(?P<pk>[0-9]+)/$', views.thread, name='thread'),
    # url(r'^(?P<pk>[0-9]+)/edit$', views.thread, name='thread'),
]