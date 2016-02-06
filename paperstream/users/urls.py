# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url, include
from . import views

urlpatterns = [
    url(r'^info/$', views.require_basic_info, name='require-basic-info'),
    url(r'^info-affiliation/$', views.require_affiliation, name='require-affiliation'),
    url(r'^email-sent/$', views.validation_sent, name='validation-sent'),
    url(r'^done/$', views.done, name='done'),
    url(r'^logout/$', views.logout, name='logout'),
    url(r'^signin/$', views.ajax_signin, name='signin'),
    url(r'^profile/$', views.profile, name='profile'),
    url(r'^settings/$', views.settings_view, name='settings'),
    url(r'^send-invite$', views.send_invite, name='send-invite'),
    url(r'^avatar/', include('avatar.urls')),
    url(r'^library/$', views.library, name='library'),
    url(r'^library/xml$', views.library_xml, name='library-xml'),
    url(r'^library/trash/$', views.library_trash, name='library-trash'),
    url(r'^library/trash/xml$', views.library_trash_xml, name='library-trash-xml'),
    url(r'^library/pins/$', views.library_likes, name='library-pins'),
    url(r'^library/pins/xml$', views.library_likes_xml, name='library-pins-xml'),
    url(r'^update-name$', views.update_name, name='update-name'),
    url(r'^update-title$', views.update_title, name='update-title'),
    url(r'^update-position$', views.update_position, name='update-position'),
    url(r'^update-affiliation$', views.update_affiliation, name='update-affiliation'),
    url(r'^delete$', views.delete_user, name='delete-user'),
    url(r'^update-stream-settings$', views.update_stream_settings, name='update-stream-settings'),
    url(r'^update-trend-settings$', views.update_trend_settings, name='update-trend-settings'),
    url(r'^update-email-digest-settings$', views.update_email_digest_settings, name='update-email-digest-settings'),
    url(r'^update-library$', views.update_library, name='update-library'),
    url(r'^user-init-step$', views.user_init_step, name='user-init-step'),
    url(r'^user-update-step$', views.user_update_step, name='user-update-step'),
    url(r'^paper/ban', views.ban_call, name='ban'),
    url(r'^paper/pin$', views.pin_call, name='pin'),
    url(r'^paper/add$', views.add_call, name='add'),
    url(r'^paper/trash$', views.trash_call, name='trash'),
    url(r'^paper/restore$', views.restore_call, name='restore'),
]
