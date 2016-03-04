# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging
from celery.canvas import chain
from config.celery import celery_app as app
from etalia.core.constants import NLP_TIME_LAPSE_CHOICES
from etalia.nlp.models import Model

logger = logging.getLogger(__name__)


@app.task()
def add_core(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


def embed_all_models(paper_pk):
    """Send chain task to embed paper
    """
    model_names = Model.objects\
        .filter(is_active=True)\
        .values_list('name', flat=True)
    for model_name in model_names:
        # Send task for embedding
        try:
            model_task = app.tasks['etalia.nlp.tasks.{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Model task for {model_name} not defined'.format(
                model_name=model_name))
            continue

        model_task.delay('infer_paper', paper_pk=paper_pk)


@app.task()
def failing_task():
    raise AssertionError
