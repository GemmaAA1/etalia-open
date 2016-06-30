# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db.models import F, Count
from django.contrib.auth import get_user_model

from config.celery import celery_app as app
from .models import Stream, Trend, ThreadFeed

logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def reset_all_main_streams():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        reset_stream.delay(user_pk)


@app.task()
def update_all_main_streams():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_stream.delay(user_pk)


@app.task()
def reset_all_main_trends():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        reset_trend.delay(user_pk)


@app.task()
def update_all_main_trends():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_trend.delay(user_pk)


@app.task()
def update_stream(user_pk, name='main'):
    """Async task / Update user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=name)
    # update
    feed.update()
    return user_pk


@app.task()
def reset_stream(user_pk, name='main'):
    """Async task / Reset user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=name)
    # reset
    feed.clear_all()
    # update
    feed.update()
    return user_pk


@app.task()
def update_trend(user_pk, name='main'):
    trend, _ = Trend.objects.get_or_create(user_id=user_pk, name=name)
    trend.update()
    return user_pk


@app.task()
def reset_trend(user_pk, name='main'):
    trend, _ = Trend.objects.get_or_create(user_id=user_pk, name=name)
    trend.update()
    return user_pk


@app.task()
def update_threadfeed(user_pk, name='main'):
    df, _ = ThreadFeed.objects.get_or_create(user_id=user_pk, name=name)
    df.update()
    return user_pk


@app.task()
def reset_threadfeed(user_pk, name='main'):
    df, _ = ThreadFeed.objects.get_or_create(user_id=user_pk, name=name)
    df.update()
    return user_pk


