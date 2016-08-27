# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^my-papers/$', views.my_papers, name='my_papers'),
    url(r'^trends/$', views.papers, name='papers'),
    url(r'^papers/(?P<pk>[0-9]+)/$', views.paper, name='paper'),
    url(r'^papers/(?P<slug>[a-zA-Z0-9-]+)_(?P<pk>[0-9]+)/$', views.paper_slug, name='paper-slug'),
]
