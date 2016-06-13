# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf.urls import include, patterns, url
from django.contrib import admin
from django.conf import settings

from etalia.core.api.router import api_urls
from etalia.feeds.views import my_feeds
from etalia.library.views import my_papers

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
        url(r'^feeds/', my_feeds, name='my_feeds'),
        url(r'^papers/', my_papers, name='my_papers'),
        url(r'^threads/', include('etalia.threads.urls', namespace='threads')),
        url(r'^user/', include('etalia.users.urls', namespace='user')),
        url(r'^user/', include('social.apps.django_app.urls', namespace='social')),
        url(r'^admin/', include(admin.site.urls)),
        url(r'^messages/', include('messages_extends.urls', namespace='message_extends')),
        # the API namespaces are defined at the app api.urls level
        url(r'^api/v1/', include(api_urls(), namespace='api')),
    ]


if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )