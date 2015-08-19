from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.feed_view, name='feed'),
    url(r'^create-feed$', views.create_feed_view, name='create-feed'),
    url(r'^(?P<feed_name>(?!create-new)[\w-]+)$', views.feed_view, name='feed'),
    url(r'^(?P<feed_name>(?!create-new)[\w-]+)/modify$', views.modify_feed_view,
        name='modify-feed'),
    url(r'^(?P<feed_name>(?!create-new)[\w-]+)/update$', views.update_feed_view,
        name='update-feed'),
]