from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^update-feed/(?P<pk>[0-9]+)/$', views.async_update_feed,
        name='update-feed'),
]