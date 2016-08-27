# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views


urlpatterns = [
    url(r'^my-threads/$', views.my_threads, name='my_threads'),
    url(r'^threads/$', views.threads, name='threads'),
    url(r'^threads/(?P<pk>[0-9]+)/$', views.thread, name='thread'),
    url(r'^threads/(?P<slug>[a-zA-Z0-9-]+)_(?P<pk>[0-9]+)/$', views.thread_slug, name='thread-slug'),
]