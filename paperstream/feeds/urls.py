from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.feed_view, name='feed'),
    url(r'^(?P<feed_name>[\w-]+)$', views.feed_view, name='feed'),
    url(r'^update-feed/(?P<pk>[0-9]+)/$', views.update_feed,
        name='update-feed'),
]