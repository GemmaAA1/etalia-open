# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app

from .models import Stats
from etalia.nlp.models import Model

import logging

logger = logging.getLogger(__name__)


@app.task()
def update_stats():
    """Create a new line of stats for the library"""
    stats = Stats.objects.create()
    stats.update()


def embed_paper(paper_pk):
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
