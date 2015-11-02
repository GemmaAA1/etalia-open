# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model

from config.celery import celery_app as app
from .models import Stream, Trend
logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def update_stream(user_pk, feed_name='main', restrict_journal=False):
    """Async task / Update user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=feed_name)
    if feed_name == 'main':
        # add all matches
        user = User.objects.get(pk=user_pk)
        feed.add_papers_seed(user.lib.papers.all())
    # update
    feed.update(restrict_journal=restrict_journal)

    return user_pk


@app.task()
def update_trend(user_pk):
    df, _ = Trend.objects.get_or_create(user_id=user_pk)
    df.update()

    return user_pk