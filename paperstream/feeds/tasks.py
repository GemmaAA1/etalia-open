# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model

from config.celery import celery_app as app
from .models import UserFeed, DiscoverFeed
logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def update_main(user_pk):
    """Async task / Init main default feed"""
    user = User.objects.get(pk=user_pk)

    # create/update main feed
    feed, _ = UserFeed.objects.get_or_create(user_id=user_pk, name='main')
    # add all papers
    feed.add_papers_seed(user.lib.papers.all())
    # update
    update_feed.delay(feed.pk)

    # Create/update Discover stream
    df, _ = DiscoverFeed.objects.get_or_create(user_id=user_pk)
    update_discover.delay(df.pk)

    return user_pk


@app.task()
def update_feed(pk):
    """Async task / Update user feed"""
    feed = UserFeed.objects.get(pk=pk)
    feed.update()


@app.task()
def update_discover(pk):
    df = DiscoverFeed.objects.get(pk=pk)
    df.update()




