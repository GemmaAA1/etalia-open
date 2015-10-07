# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app

from .models import Stats


@app.task()
def update_stats():
    stats = Stats.objects.first()
    stats.update()
