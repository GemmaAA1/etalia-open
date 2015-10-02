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
def embed_paper(pk, model_name):
    model = Model.objects.get(name=model_name)
    model.infer_paper(paper_pk=pk)
