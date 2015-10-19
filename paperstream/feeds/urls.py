# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.feed_main, name='main'),
    url(r'^discover$', views.trend_view, name='trend'),
    url(r'^create-feed$', views.create_feed_view, name='create-feed'),
    url(r'^(?P<name>[\w-]+)/$', views.feed_view, name='feed'),
    url(r'^(?P<name>[\w-]+)/modify$', views.modify_feed_view,
        name='modify-feed'),
    url(r'^(?P<name>[\w-]+)/delete$', views.delete_feed_view,
        name='delete-feed'),
    url(r'^(?P<name>[\w-]+)/update$', views.update_feed_view,
        name='update-feed'),
    url(r'^(?P<name>[\w-]+)/user-feed-message$', views.ajax_user_feed_message,
        name='user-feed-message'),
]
