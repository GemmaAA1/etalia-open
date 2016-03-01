# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf.urls import include, patterns, url
from django.contrib import admin
from django.conf import settings

if settings.INVITE_MODE:
    urlpatterns = [
        url(r'^', include('etalia.invite.urls', namespace='invite')),
        url(r'^library/', include('etalia.library.urls', namespace='library')),
    ]
else:
    urlpatterns = [
        url(r'^', include('etalia.core.urls', namespace='core')),
        url(r'^library/', include('etalia.library.urls', namespace='library')),
        url(r'^feed/', include('etalia.feeds.urls', namespace='feeds')),
        url(r'^user/', include('etalia.users.urls', namespace='user')),
        url(r'^user/', include('social.apps.django_app.urls', namespace='social')),
        url(r'^admin/', include(admin.site.urls)),
        url(r'^messages/', include('messages_extends.urls', namespace='message_extends')),
    ]


if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )