# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import glob
import logging
from celery import chain, Task

from config.celery import celery_app as app

from django.conf import settings

from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES
from paperstream.core.utils import db_table_exists

from .models import Model, MostSimilar

logger = logging.getLogger(__name__)


class EmbedPaperTask(Task):
    """Abstract Embedding Paper task for model
    Use to load model in __init__ so that it is cached for subsequent call
    to task
    """
    ignore_result = False
    model_name = None
    _model = None

    def __init__(self, *args, **kwargs):
        if 'model_name' in kwargs:
            model_name = kwargs['model_name']
            # check in model_name is known
            choices = [model['name'] for model in
                       Model.objects.all().values('name')]
            if model_name in choices:
                self.model_name = model_name
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))

    @property
    def model(self):
        if self._model is None:
            self._model = Model.objects.load(name=self.model_name)
            return self._model
        # if Model has been modified and currently not uploading, reload
        model_now = Model.objects.get(name=self.model_name)
        last_modified = model_now.modified
        upload_state = model_now.upload_state
        if not self._model.modified == last_modified and upload_state == 'IDL':
            self._model = Model.objects.load(name=self.model_name)

        return self._model

    def run(self, paper_pk, **kwargs):
        return self.model.infer_paper(paper_pk)


class MostSimilarTask(Task):
    """Abstract task for MostSimilar model
    Use to load MostSimilar instance in __init__ so that it is cached for subsequent call
    to task
    """
    ignore_result = False
    model_name = None
    _ms = None

    def __init__(self, *args, **kwargs):
        if 'model_name' in kwargs:
            model_name = kwargs['model_name']

            # check if model_name is known
            choices = [model['name'] for model in
                       Model.objects.all().values('name')]
            if model_name in choices:
                self.model_name = model_name
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))

    @property
    def ms(self):
        if self._ms is None:
            self._ms = MostSimilar.objects.load(model__name=self.model_name)
            return self._ms
        # if MostSimilar has been modified and MostSimilar not uploading, reload
        ms_now = MostSimilar.objects.get(model__name=self.model_name)
        last_modified = ms_now.modified
        upload_state = ms_now.upload_state
        if not self._ms.modified == last_modified and upload_state == 'IDL':
            # remove local
            rm_files = glob.glob(
                os.path.join(settings.NLP_MS_PATH, '{name}.ms*'.format(
                    name=self._ms.name)))
            for file in rm_files:
                os.remove(file)

            self._ms = MostSimilar.objects.load(model__name=self.model_name)
        return self._ms

    def run(self, *args, **kwargs):
        return self.ms.tasks(*args, **kwargs)

# Instantiate all task for model and Mo
def register_nlp_tasks():
    # Create embedding task from model
    model_names = Model.objects.all().values_list('name', flat=True)

    for model_name in model_names:
        cls = EmbedPaperTask(model_name=model_name)
        app.task(cls, name='paperstream.nlp.tasks.{model_name}'.format(
            model_name=model_name))

        cls = MostSimilarTask(model_name=model_name)
        app.task(cls, name='paperstream.nlp.tasks.mostsimilar_{model_name}'.format(
            model_name=model_name))


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

# If table exist register task
# NB: this avoid failing when running manage.py migrate from scratch
if db_table_exists('nlp_model'):
    register_nlp_tasks()
