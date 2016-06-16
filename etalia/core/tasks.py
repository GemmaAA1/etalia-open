# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from config.celery import celery_app as app

logger = logging.getLogger(__name__)


@app.task()
def add_core(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


@app.task()
def failing_task():
    """For email error testing"""
    raise AssertionError