# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging
from celery.canvas import chain
from config.celery import celery_app as app
from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES
from paperstream.nlp.models import Model

logger = logging.getLogger(__name__)


@app.task()
def add_core(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


def embed_all_models(paper_pk):
    """Send chain task to embed paper
    """
    model_names = Model.objects.all().values_list('name', flat=True)
    for model_name in model_names:
        # Send task for embedding
        try:
            embed_task = app.tasks['paperstream.nlp.tasks.{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Embeding task for {model_name} not defined'.format(
                model_name=model_name))
            continue

        embed_task.apply_async(args=(paper_pk, ))
