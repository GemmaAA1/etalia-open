# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from .models import Model, PaperEngine, ThreadEngine
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

# NOTE:
# tasks related to Model or Engines are registered in celery.py because they
# are host dependant


@app.task()
def paperengine_update_all():
    models = Model.objects.filter(is_active=True)
    for model in models:
        pe = PaperEngine.objects.load(model=model,
                                      is_active=True)
        pe.update()


@app.task()
def paperengine_full_update_all():
    models = Model.objects.filter(is_active=True)
    for model in models:
        pe = PaperEngine.objects.load(model=model,
                                      is_active=True)
        pe.full_update()
        pe.activate()

@app.task()
def threadengine_update_all():
    models = Model.objects.filter(is_active=True)
    for model in models:
        te = ThreadEngine.objects.load(model=model,
                                       is_active=True)
        te.update()

@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


