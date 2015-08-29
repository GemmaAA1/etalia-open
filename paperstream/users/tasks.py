# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model
from celery.canvas import chain

from config.celery import celery_app as app
from paperstream.feeds.tasks import init_main_feed

logger = logging.getLogger(__name__)

User = get_user_model()


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


def init_user(user_pk, provider_name):
    """Task init user / Chain user library update, and main feed initialization
    """
    task = chain(update_lib.s(user_pk, provider_name),
                 init_main_feed.s())

    task.delay()
    return user_pk
