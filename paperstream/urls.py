from django.conf.urls import include, url
from django.contrib import admin
from base import views as base_views

urlpatterns = [
    # Examples:
    url(r'^$', base_views.home, name='home'),
    url(r'^library/', include('library.urls')),
    # url(r'^accounts/', include('accounts.urls')),
    # url(r'^accounts/', include('allauth.urls')),
    url(r'^admin/', include(admin.site.urls)),
]
