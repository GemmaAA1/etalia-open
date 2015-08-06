import logging

from celery import chain
from config.celery import celery_app as app

from feeds.models import UserFeed

logger = logging.getLogger(__name__)

@app.task()
def update_feed(feed_pk):

    # get feed
    feed = UserFeed.objects.get(pk=feed_pk)

    # update
    feed.update()

    return feed_pk