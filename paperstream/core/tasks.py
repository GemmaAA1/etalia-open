# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging
from celery.canvas import chain
from config.celery import celery_app as app
from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES
from paperstream.nlp.models import Model

logger = logging.getLogger(__name__)


@app.task
def add(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


def embed_all_models_and_find_neighbors(paper_pk):
    """Send chain task to embed paper and populate neighbors
    """
    model_names = Model.objects.all().values_list('name', flat=True)
    for model_name in model_names:
        # Send task for embedding
        try:
            embed_task = app.tasks['paperstream.nlp.tasks.embed_paper_{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Embeding task for {model_name} not defined'.format(
                model_name=model_name))
            continue

        # Send task for time_lapse related LSHs
        for time_lapse, _ in NLP_TIME_LAPSE_CHOICES:
            try:
                lsh_task = app.tasks['paperstream.nlp.tasks.lsh_{model_name}_{time_lapse}'
                    .format(model_name=model_name, time_lapse=time_lapse)]
            except KeyError:
                logger.error('LSH task for {model_name}/{time_lapse} not defined'
                             .format(model_name=model_name,
                                     time_lapse=time_lapse))
                continue

            task = chain(embed_task.s(paper_pk),
                         lsh_task.s('populate_neighbors'))
            task.apply_async()
