# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model

from config.celery import celery_app as app
from .models import Stream, Trend
logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def reset_all_main_streams():
    us_pk = User.objects.all().values_list('pk', flat=True)
    for user_pk in us_pk:
        reset_stream.delay(user_pk)


@app.task()
def update_all_main_streams():
    us_pk = User.objects.all().values_list('pk', flat=True)
    for user_pk in us_pk:
        update_stream.delay(user_pk)


@app.task()
def reset_all_main_trends():
    us_pk = User.objects.all().values_list('pk', flat=True)
    for user_pk in us_pk:
        reset_trend.delay(user_pk)


@app.task()
def update_all_main_trends():
    us_pk = User.objects.all().values_list('pk', flat=True)
    for user_pk in us_pk:
        update_trend.delay(user_pk)


@app.task()
def update_stream(user_pk, stream_name='main', restrict_journal=False):
    """Async task / Update user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=stream_name)
    user = User.objects.get(pk=user_pk)
    if stream_name == 'main':
        # add all seeds
        feed.add_papers_seed(user.lib.papers.all())
    # update
    feed.update(restrict_journal=restrict_journal)
    return user_pk


@app.task()
def reset_stream(user_pk, stream_name='main', restrict_journal=False):
    """Async task / Reset user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=stream_name)
    user = User.objects.get(pk=user_pk)
    if stream_name == 'main':
        # add all seeds
        feed.add_papers_seed(user.lib.papers.all())
    # reset
    feed.clear_all()
    if hasattr(user, 'streamlayout'):
        user.streamlayout.delete()
    # update
    feed.update(restrict_journal=restrict_journal)
    return user_pk


@app.task()
def update_trend(user_pk, trend_name='main'):
    df, _ = Trend.objects.get_or_create(user_id=user_pk, name=trend_name)
    df.update()
    return user_pk


@app.task()
def reset_trend(user_pk, trend_name='main'):
    df, _ = Trend.objects.get_or_create(user_id=user_pk, name=trend_name)
    user = User.objects.get(pk=user_pk)
    # reset filter
    if hasattr(user, 'trendlayout'):
        user.trendlayout.delete()
    # update
    df.update()
    return user_pk