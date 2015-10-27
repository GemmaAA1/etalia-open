# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf.urls import include, url
from django.contrib import admin

urlpatterns = [
    # Examples:
    url(r'^', include('paperstream.core.urls', namespace='core')),
    url(r'^library/', include('paperstream.library.urls', namespace='library')),
    url(r'^feed/', include('paperstream.feeds.urls', namespace='feeds')),
    url(r'^user/', include('paperstream.users.urls', namespace='user')),
    url(r'^user/', include('social.apps.django_app.urls', namespace='social')),
    url(r'^admin/', include(admin.site.urls)),
]
