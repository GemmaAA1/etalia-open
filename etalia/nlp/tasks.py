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
    pe_ids = PaperEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for peid in pe_ids:
        pe = PaperEngine.objects.load(id=peid, is_active=True)
        pe.update()


@app.task()
def paperengine_full_update_all():
    pe_ids = PaperEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for peid in pe_ids:
        pe = PaperEngine.objects.load(id=peid, is_active=True)
        pe.full_update()


@app.task()
def threadengine_update_all():
    te_ids = ThreadEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for teid in te_ids:
        te = ThreadEngine.objects.load(id=teid, is_active=True)
        te.update()


@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


