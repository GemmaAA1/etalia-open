# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from .models import Model, MostSimilar
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

# PLEASE NOTE:
# tasks related to embeding of paper and mostsimilar model are registered
# in celery.py because they are host dependant

@app.task()
def mostsimilar_update_all():
    models = Model.objects.filter(is_active=True)
    for model in models:
        ms = MostSimilar.objects.load(model=model,
                                      is_active=True)
        ms.update()

@app.task()
def mostsimilar_full_update_all():
    models = Model.objects.filter(is_active=True)
    for model in models:
        ms = MostSimilar.objects.load(model=model,
                                      is_active=True)
        ms.full_update()
        ms.activate()

@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


