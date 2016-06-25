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
    url(r'^update-name$', views.update_name, name='update-name'),
    url(r'^update-title$', views.update_title, name='update-title'),
    url(r'^update-position$', views.update_position, name='update-position'),
    url(r'^update-affiliation$', views.update_affiliation, name='update-affiliation'),
    url(r'^delete$', views.delete_user, name='delete-user'),
    url(r'^update-fingerprint-settings$', views.update_fingerprint_settings, name='update-fingerprint-settings'),
    url(r'^update-stream-settings$', views.update_stream_settings, name='update-stream-settings'),
    url(r'^update-trend-settings$', views.update_trend_settings, name='update-trend-settings'),
    url(r'^update-email-digest-settings$', views.update_email_digest_settings, name='update-email-digest-settings'),
    url(r'^update-library$', views.update_library, name='update-library'),
    url(r'^user-init-check$', views.user_init_check, name='user-init-check'),
    url(r'^user-update-stream-check$', views.user_update_stream_check, name='user-update-stream-check'),
    url(r'^user-update-trend-check$', views.user_update_trend_check, name='user-update-trend-check'),
    url(r'^user-update-settings-check$', views.user_update_settings_check, name='user-update-settings-check'),
]
