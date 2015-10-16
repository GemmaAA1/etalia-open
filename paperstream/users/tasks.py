# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model
from django.utils import timezone
from celery.canvas import chain

from config.celery import celery_app as app
from paperstream.feeds.tasks import update_feed, update_discover
from paperstream.library.models import Paper

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


@app.task()
def add_paper_to_lib(user_pk, provider_name, paper_pk):
    """Async add paper to usr library"""
    # get user
    user = User.objects.get(pk=user_pk)

    # get social
    social = user.social_auth.get(provider=provider_name)

    # get backend
    backend = social.get_backend_instance()

    # build session
    session = backend.get_session(social, user)

    # Get paper
    paper = Paper.objects.get(pk=paper_pk)

    # push paper to lib
    err = backend.add_paper(session, paper)

    # add paper locally
    backend.associate_paper(paper, user, {'created': timezone.now().date()})
    backend.associate_journal(paper.journal, user)

    return err


def init_user(user_pk, provider_name):
    """Task init user / Chain user library update, and feed initialization
    """
    task = chain(update_lib.s(user_pk, provider_name),
                 update_feed.s('main'),
                 update_discover.s())

    # task.delay()
    task()
    return user_pk
