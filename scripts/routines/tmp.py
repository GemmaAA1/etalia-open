# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from etalia.feeds.tasks import reset_stream, reset_trend
from django.contrib.auth import get_user_model


User = get_user_model()
us_pk = User.objects.all().values_list('pk', flat=True)
user_pk = us_pk[0]

reset_trend(user_pk)


tmp = {}


def add_some(tmp):
    # tmp['a'] = 1
    tmp = tmp + 1


def filter_tmp(qs):
    qs.