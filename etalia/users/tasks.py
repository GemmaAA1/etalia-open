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
        provider_names = user.social_auth.all().values_list('provider', flat=True)
        for provider_name in provider_names:
            update_lib.delay(user.id, provider_name)


@app.task()
def update_lib(user_pk, provider_name):
    """Async task for updating user library"""
    # get user
    user = User.objects.get(pk=user_pk)

    # get social
    social = user.social_auth.get(provider=provider_name)

    # get backend
    backend = social.get_backend_instance()

    # build session
    session = backend.get_session(social, user)

    # update lib
    backend.update_lib(user, session)

    return user_pk


@app.task()
def init_step(user_pk, step):
    """Set flag is_init to True"""
    user = User.objects.get(pk=user_pk)
    user.init_step = step
    user.save()

    return user_pk


def init_user(user_pk, provider_name):
    """Task init user / Chain user library update, and feed initialization
    """
    task = chain(
        init_step.s(user_pk, 'LIB'),
        update_lib.s(provider_name),
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
