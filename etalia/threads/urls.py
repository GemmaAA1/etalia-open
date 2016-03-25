# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views


urlpatterns = [
    url(r'^$', views.my_threads, name='my_threads'),
    url(r'^new$', views.new_thread, name='new_thread'),
    url(r'^(?P<pk>[0-9]+)/$', views.thread, name='thread'),
    url(r'^(?P<pk>[0-9]+)/edit', views.update_thread, name='update_thread'),
    url(r'^(?P<pk>[0-9]+)/posts/new$', views.new_post, name='new_post'),
    url(r'^posts/(?P<pk>[0-9]+)/edit', views.edit_post, name='edit_post'),
    url(r'^posts/(?P<pk>[0-9]+)/delete', views.delete_post, name='delete_post'),
    url(r'^posts/(?P<pk>[0-9]+)/comments/new$', views.new_comment, name='new_comment'),
    url(r'^comments/(?P<pk>[0-9]+)/edit', views.edit_comment, name='edit_comment'),
    url(r'^comments/(?P<pk>[0-9]+)/delete$', views.delete_comment, name='delete_comment'),
    url(r'^state/', views.thread_state, name='state'),
]