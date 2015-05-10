from django.conf.urls import include, url
from django.contrib import admin
from base import views

urlpatterns = [
    # Examples:
    url(r'^$', include('base.urls', namespace='base')),
    url(r'^library/', include('library.urls', namespace='library')),
    # url(r'^accounts/', include('accounts.urls')),
    # url(r'^accounts/', include('allauth.urls')),
    url(r'^admin/', include(admin.site.urls)),
]
