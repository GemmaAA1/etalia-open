# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf.urls import include, patterns, url
from django.contrib import admin
from django.views.generic import TemplateView
from django.conf import settings

from etalia.core.api.router import api_urls


urlpatterns = [
    url(r'^', include('etalia.core.urls', namespace='core')),
    url(r'^', include('etalia.library.urls', namespace='library')),
    url(r'^', include('etalia.feeds.urls', namespace='feeds')),
    url(r'^', include('etalia.threads.urls', namespace='threads')),
    url(r'^user/', include('etalia.users.urls', namespace='user')),
    url(r'^user/', include('social.apps.django_app.urls', namespace='social')),
    url(r'^admin/', include(admin.site.urls)),
    # the API namespaces are defined at the app api.urls level
    url(r'^api/v1/', include(api_urls(), namespace='api')),
    url(r'^robots\.txt$', TemplateView.as_view(
        template_name="robots.txt",
        content_type={'mimetype': 'text/plain'})),
]


if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )