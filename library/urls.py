from django.conf.urls import url
from library import views

urlpatterns = [
    url(r'^$', views.view_library, name='view_library'),
    url(r'^journals/$', views.view_journals, name='view_journals'),
    url(r'^journal/(?P<id>[0-9]+)/$', views.view_journal, name='view_journal'),
    url(r'^papers/$', views.view_papers, name='view_papers'),
    url(r'^paper/(?P<id>[0-9]+)/$', views.view_paper, name='view_paper'),
]
