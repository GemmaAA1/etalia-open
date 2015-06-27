from django.conf.urls import url, include
from . import views

urlpatterns = [
    url(r'^info/$', views.require_basic_info, name='require-basic-info'),
    url(r'^info-affiliation/$', views.require_affiliation,
        name='require-affiliation'),
    url(r'^email-sent/$', views.validation_sent, name='validation-sent'),
    url(r'^done/$', views.done, name='done'),
    url(r'^logout/$',views.logout, name='logout'),
    # url(r'^signin/$',views.signin, name='signin'),
    url(r'^signin/$',views.ajax_signin, name='signin'),
    url(r'^library/$',views.library, name='library'),
    url(r'^update-basic-info/$',views.ajax_update_basic_info,
        name='update-basic-info'),
    url(r'^update-affiliation/$',views.ajax_update_affiliation,
        name='update-affiliation'),
    url(r'^user-lib-count-papers/$', views.ajax_user_lib_count_papers,
        name='user-lib-count-papers'),
    url(r'^update-user-lib/$', views.async_update_user_lib,
        name='update-user-lib'),
]
