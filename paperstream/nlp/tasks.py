# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from .models import Model, MostSimilar
from .tasks_class import EmbedPaperTask, MostSimilarTask
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

# NB: tasks related to embeding of paper and mostsimilar are register in celery.py

@app.task()
def mostsimilar_update_all():
    models = Model.objects.all()
    for model in models:
        ms = MostSimilar.objecs.load(model=model)
        ms.update()


@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y


@app.task()
def embed_papers_batch(pks, model_name):
    pks = list(pks)
    embed_task = app.tasks['paperstream.nlp.tasks.{model_name}'.format(
        model_name=model_name)]
    for pk in pks[:-1]:
        embed_task.apply_async(args=(pk, ))
    # we wait for the results of the last one
    pk = pks[-1]
    embed_task = app.tasks['paperstream.nlp.tasks.{model_name}'.format(
        model_name=model_name)]
    res = embed_task.apply_async(args=(pk, ))
    return res.get()


def embed_papers(pks, model_name, batch_size=10000):
    pks = list(pks)
    # split in batchs
    pks_batch = [pks[i*batch_size:(i+1)*batch_size] for i in range(len(pks)//batch_size)]
    pks_batch.append(pks[len(pks)//batch_size*batch_size:])

    for batch in pks_batch:
        embed_papers_batch.delay(batch, model_name)

