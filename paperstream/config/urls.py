from django.conf.urls import include, url
from django.contrib import admin

urlpatterns = [
    # Examples:
    url(r'^$', include('core.urls', namespace='core')),
    url(r'^library/', include('library.urls', namespace='library')),
    url(r'^users/', include('users.urls', namespace='users')),
    url(r'^users/', include('social.apps.django_app.urls', namespace='social')),
    url(r'^admin/', include(admin.site.urls)),
]
