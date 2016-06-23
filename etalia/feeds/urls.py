# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    # url(r'^$', views.stream_main, name='main'),
    url(r'^$', views.my_feeds, name='my_feeds'),
    url(r'stream/(?P<name>[\w-]+)/update$', views.update_stream_view, name='update-stream'),
    url(r'stream/(?P<name>[\w-]+)/reset$', views.reset_stream_view, name='reset-stream'),
    url(r'trend/(?P<name>[\w-]+)/update$', views.update_trend_view, name='update-trend'),
    url(r'trend/(?P<name>[\w-]+)/reset$', views.reset_trend_view, name='reset-trend'),
]
