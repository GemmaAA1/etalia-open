from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.feed_view, name='main'),
    url(r'^create-feed$', views.create_feed_view, name='create-feed'),
    url(r'^(?P<pk>[\d]+)/$', views.feed_view, name='feed'),
    url(r'^(?P<pk>[\d]+)/modify$', views.modify_feed_view,
        name='modify-feed'),
    url(r'^(?P<pk>[\d]+)/delete$', views.delete_feed_view,
        name='delete-feed'),
    url(r'^(?P<pk>[\d]+)/update$', views.update_feed_view,
        name='update-feed'),
    url(r'^(?P<pk>[\d]+)/user-feed-message$', views.ajax_user_feed_message,
        name='user-feed-message'),
]
