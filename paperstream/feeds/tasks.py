import logging

from django.contrib.auth import get_user_model

from config.celery import celery_app as app
from .models import UserFeed
logger = logging.getLogger(__name__)

User = get_user_model()

@app.task()
def init_main_feed(user_pk):
    """Async task / Init main default feed"""
    user = User.objects.get(pk=user_pk)

    # get default ('main') feed
    feed = UserFeed.objects.get(user_id=user_pk, name='main')

    # add all papers
    feed.add_papers_seed(user.lib.papers.all())

    # update
    feed.update()

    return user_pk


@app.task()
def update_feed(pk):
    """Async task / Update user feed"""
    feed = UserFeed.objects.get(pk=pk)
    feed.update()




