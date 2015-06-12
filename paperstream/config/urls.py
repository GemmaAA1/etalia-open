from django.conf.urls import include, url
from django.contrib import admin

from core.views import test

urlpatterns = [
    # Examples:
    url(r'^$', include('core.urls', namespace='core')),
    url(r'^test/$', test, name='test'),
    url(r'^library/', include('library.urls', namespace='library')),
    url(r'^feed/', include('feeds.urls', namespace='feeds')),
    url(r'^user/', include('users.urls', namespace='user')),
    url(r'^user/', include('social.apps.django_app.urls', namespace='social')),
    url(r'^admin/', include(admin.site.urls)),
]
