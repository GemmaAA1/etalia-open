# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import glob
from celery import Task

from django.db.models.query import QuerySet
from django.conf import settings

from .models import Model, PaperEngine, ThreadEngine


class EmbedPaperTask(Task):
    """Abstract Embedding Paper task for model
    Use to load model in __init__ so that it is cached for subsequent call
    to task
    """

    abstract = True
    ignore_result = False
    model_name = None
    _model = None

    def __init__(self, *args, **kwargs):
        if 'model_name' in kwargs:
            model_name = kwargs.get('model_name')
            # check in model_name is known
            choices = [model['name'] for model in
                       Model.objects.all().values('name')]
            if model_name in choices:
                self.model_name = model_name
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))

    def load(self):
        self._model = Model.objects.load(name=self.model_name)
        return self._model

    @property
    def model(self):
        if self._model is None:
            return self.load()
        # if Model has been modified and currently not uploading, reload
        model_now = Model.objects.get(name=self.model_name)
        last_modified = model_now.modified
        upload_state = model_now.upload_state
        if not self._model.modified == last_modified and upload_state == 'IDL':
            return self.load()

        return self._model


class EngineTask(Task):
    """Abstract task for Engine model

    Use to load Engine (PaperEngine, ThreadEnging) type of instance in __init__
    so that it is cached for subsequent call to task (avoid overhead of loading
    data for each task)
    """

    abstract = True
    engine_class = None
    engine_id = None
    ignore_result = False

    _engine = None

    def __init__(self, *args, **kwargs):
        if 'engine_id' in kwargs:
            # check in engine_id is known
            engine_id = kwargs.get('engine_id')
            choices = self.engine_class.objects.filter(is_active=True).values_list('id', flat=True)
            if engine_id in choices:
                self.engine_id = engine_id
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))

    def load(self):
        self._engine = self.engine_class.objects.load(
            id=self.engine_id,
            is_active=True)
        return self._engine

    def get(self):
        return self.engine_class.objects.get(id=self.engine_id, is_active=True)

    @property
    def engine(self):
        if self._engine is None:
            return self.load()
        # if Engine has been modified and Engine not uploading/downloading, or
        # has been deactivated => reload
        engine_now = self.get()
        last_modified = engine_now.modified
        upload_state = engine_now.upload_state
        active_now = engine_now.is_active
        if (not active_now) or \
                (upload_state == 'IDL' and not self._engine.modified == last_modified):
            # remove local
            rm_files = glob.glob(
                os.path.join(self._engine.PATH,
                             '{0}.*'.format(self._engine.name)))
            for file in rm_files:
                os.remove(file)
            return self.load()

        return self._engine


class PaperEngineTask(EngineTask):
    abstract = True
    engine_class = PaperEngine


class ThreadEngineTask(EngineTask):
    abstract = True
    engine_class = ThreadEngine
