# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.library, name='library'),
    url(r'^journals/$', views.journals, name='journals'),
    url(r'^journal/(?P<pk>[0-9]+)/$', views.journal, name='journal'),
    url(r'^journal/(?P<slug>[a-zA-Z0-9-]+)-(?P<pk>[0-9]+)/$', views.journal_slug, name='journal-slug'),
    url(r'^paper/(?P<pk>[0-9]+)/$', views.paper, name='paper'),
    url(r'^paper/(?P<pk>[0-9]+)/time=(?P<time_lapse>[\w]+)$', views.paper_time, name='paper-time'),
    url(r'^paper/(?P<slug>[a-zA-Z0-9-]+)-(?P<pk>[0-9]+)/$', views.paper_slug, name='paper-slug'),
    url(r'^paper/(?P<slug>[a-zA-Z0-9-]+)-(?P<pk>[0-9]+)/time=(?P<time_lapse>[\w]+)$', views.paper_slug, name='paper-slug-time'),
]
