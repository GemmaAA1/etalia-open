# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import glob
import logging
from celery import Task

from django.conf import settings

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
    abstract = True

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
    abstract = True

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
            # remove local
            rm_files = glob.glob(
                os.path.join(settings.NLP_MS_PATH, '{name}.ms*'.format(
                    name=self.model_name)))
            for file in rm_files:
                os.remove(file)
            self._ms = MostSimilar.objects.load(model__name=self.model_name)
            return self._ms
        # if MostSimilar has been modified and MostSimilar not uploading, reload
        ms_now = MostSimilar.objects.get(model__name=self.model_name)
        last_modified = ms_now.modified
        upload_state = ms_now.upload_state
        if upload_state == 'IDL' and not self._ms.modified == last_modified:
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
