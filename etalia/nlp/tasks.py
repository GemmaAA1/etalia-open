# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from .models import Model, MostSimilar
from .tasks_class import EmbedPaperTask, MostSimilarTask
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
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


def embed_papers(pks, model_name, batch_size=1000):
    try:
        model_task = app.tasks['etalia.nlp.tasks.{model_name}'.format(
            model_name=model_name)]
    except KeyError:
        logger.error('Embeding task for {model_name} not defined'.format(
            model_name=model_name))
        raise KeyError
    pks = list(pks)
    nb_papers = len(pks)
    nb_batches = nb_papers // batch_size
    pks_batched = [pks[i*batch_size:(1+i)*batch_size] for i in range(nb_batches)]
    pks_batched.append(pks[nb_batches * batch_size:])

    for batch in pks_batched:
        model_task.delay('infer_papers', batch)