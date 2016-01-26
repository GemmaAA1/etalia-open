# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    # url(r'^$', views.stream_main, name='main'),
    url(r'stream/$', views.stream_view, name='home'),
    url(r'stream/$', views.stream_view, name='stream'),
    url(r'stream/filter$', views.stream_view2, name='stream_filter'),
    url(r'stream/(?P<name>[\w-]+)/update$', views.update_stream_view,
        name='update-stream'),
    url(r'stream/(?P<name>[\w-]+)/reset$', views.reset_stream_view,
        name='reset-stream'),
    url(r'trend/$', views.trend_view, name='trend'),
    url(r'trend/filter$', views.trend_view2, name='trend_filter'),
    url(r'trend/(?P<name>[\w-]+)/update$', views.update_trend_view,
        name='update-trend'),
    url(r'trend/(?P<name>[\w-]+)/reset$', views.reset_trend_view,
        name='reset-trend'),
    # url(r'^create-feed$', views.create_feed_view, name='create-feed'),
    # url(r'^(?P<name>[\w-]+)/$', views.stream_view, name='stream'),
    # url(r'^(?P<name>[\w-]+)/modify$', views.modify_feed_view,
    #     name='modify-feed'),
    # url(r'^(?P<name>[\w-]+)/delete$', views.delete_feed_view,
    #     name='delete-feed'),
    # url(r'^(?P<name>[\w-]+)/update$', views.update_feed_view,
    #     name='update-feed'),
    # url(r'^(?P<name>[\w-]+)/user-feed-message$', views.ajax_user_feed_message,
    #     name='user-feed-message'),
]
