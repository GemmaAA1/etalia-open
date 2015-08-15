import logging

from config.celery import celery_app as app
from .models import UserFeed

logger = logging.getLogger(__name__)

@app.task()
def update_feed(pk):
    """Update user feed"""
    feed = UserFeed.objects.get(pk=pk)
    feed.update()
