# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model
from celery.canvas import chain

from config.celery import celery_app as app
from etalia.feeds.tasks import update_stream, update_trend, update_threadfeed

logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def userlib_update_all():
    us = User.objects.all()
    for user in us:
        update_lib.delay(user.id)


@app.task()
def update_lib(user_id):
    """Async task for updating user library"""
    user = User.objects.get(pk=user_id)
    user.lib.update()
    return user_id


@app.task()
def init_step(user_pk, step):
    """Set flag is_init to True"""
    user = User.objects.get(pk=user_pk)
    user.init_step = step
    user.save()

    return user_pk


def init_user(user_pk):
    """Task init user / Chain user library update, and feed initialization
    """
    task = chain(
        init_step.s(user_pk, 'LIB'),
        update_lib.s(),
        init_step.s('STR'),
        update_stream.s(),
        init_step.s('TRE'),
        update_trend.s(),
        init_step.s('THR'),
        update_threadfeed.s(),
        init_step.s('IDL'),
    )

    task()
    return user_pk
