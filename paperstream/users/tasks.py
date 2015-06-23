import logging
from django.conf import settings
from django.contrib.auth import get_user_model

from celery import chain
from config.celery import celery_app as app

from feeds.models import UserFeed

logger = logging.getLogger(__name__)

User = get_user_model()


@app.task(name='tasks.init_user')
def init_user(user_pk, provider_name):
    # chain user library update, and default feed initialization
    chain(update_lib.s(user_pk, provider_name), init_default_feed.s())()

    return user_pk

@app.task(name='tasks.update_lib')
def update_lib(user_pk, provider_name):

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


# The following task may be moved to feeds.tasks

@app.task(name='tasks.init_default_feed')
def init_default_feed(user_pk):

    # get user
    user = User.objects.get(pk=user_pk)

    # create default ('main') feed
    user_feed = UserFeed.objects.init_default_userfeed(user)

    # update
    user_feed.update()

    return user_pk

@app.task(name='tasks.update_feed')
def update_feed(user_pk, feed_name):

    # get user
    user = User.objects.get(pk=user_pk)

    # get user_feed
    user_feed = user.feed.get(name='feed_name')

    # update
    user_feed.update()

    return user_pk

