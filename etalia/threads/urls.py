# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views


urlpatterns = [
    url(r'^$', views.library, name='library'),

]