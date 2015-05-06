from django.conf.urls import url
from library import views

urlpatterns = [
    url(r'^$', views.view_library, name='view_library'),
    url(r'^papers/(\d+)/$', views.view_papers, name='view_papers'),
    url(r'^journals/(\d+)/$', views.view_journals, name='view_journals'),
]
