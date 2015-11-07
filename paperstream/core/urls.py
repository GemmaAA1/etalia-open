# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^about/$', views.about, name='about'),
    url(r'^terms/$', views.terms, name='terms'),
    url(r'^support/$', views.support, name='support'),
    url(r'^help/$', views.help, name='help'),
]
