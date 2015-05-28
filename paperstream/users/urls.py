from django.conf.urls import url, include
from . import views

urlpatterns = [
    url(r'^info/$', views.require_primary, name='require_primary'),
    url(r'^info-affiliation/$', views.require_affiliation, name='require_affiliation'),
    url(r'^email-sent/$', views.validation_sent, name='validation_sent'),
    url(r'^done/$', views.done, name='done'),
    url(r'^logout/$',views.logout, name='logout'),
]
