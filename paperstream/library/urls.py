# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.library, name='library'),
    url(r'^journals/$', views.journals, name='journals'),
    url(r'^journal/(?P<pk>[0-9]+)/$', views.journal, name='journal'),
    url(r'^paper/(?P<pk>[0-9]+)/$', views.paper, name='paper'),
    url(r'^papers/$', views.papers, name='papers'),
]
