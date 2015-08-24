import logging
from celery import chain
from config.celery import celery_app as app

from feeds.tasks import init_main_feed
from users.tasks import update_lib

logger = logging.getLogger(__name__)

@app.task()
def init_user(user_pk, provider_name):
    """Task init user / Chain user library update, and main feed initialization
    """
    chain(update_lib.s(user_pk, provider_name),
          init_main_feed.s())()

    return user_pk
