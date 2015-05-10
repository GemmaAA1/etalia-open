from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.library, name='library'),
    url(r'^journals/$', views.journals, name='journals'),
    url(r'^journal/(?P<id>[0-9]+)/$', views.journal, name='journal'),
    url(r'^papers/$', views.papers, name='papers'),
    url(r'^paper/(?P<id>[0-9]+)/$', views.paper, name='paper'),
]
