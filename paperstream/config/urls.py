from django.conf.urls import include, url
from django.contrib import admin

urlpatterns = [
    # Examples:
    url(r'^$', include('core.urls', namespace='core')),
    url(r'^library/', include('library.urls', namespace='library')),
    # url(r'^users/', include('users.urls')),
    # url(r'^accounts/', include('allauth.urls')),
    url(r'^admin/', include(admin.site.urls)),
]