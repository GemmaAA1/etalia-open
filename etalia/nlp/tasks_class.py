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
    ignore_result = False
    model_name = None
    _model = None
    init = False

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

        # init task
        self.init = kwargs.get('init', False)
        if self.init:
            self._model = Model.objects.load(name=self.model_name)

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

    def run(self, *args, **kwargs):
        return self.model.tasks(*args, **kwargs)


class EngineTask(Task):
    """Abstract task for Engine model

    Use to load Engine (PaperEngine, ThreadEnging) type of instance in __init__
    so that it is cached for subsequent call to task (avoid overhead of loading
    data for each task)
    """
    engine_class = None
    engine_id = None
    ignore_result = False
    init = False

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

        # init task
        self.init = kwargs.get('init', False)
        if self.init:
            self._engine = self.engine_class.objects.load(
                id=self.engine_id,
                is_active=True)

    @property
    def engine(self):
        if self._engine is None:
            self._engine = self.engine_class.objects.load(
                id=self.engine_id,
                is_active=True)
            return self._engine
        # if Engine has been modified and Engine not uploading/downloading, or
        # has been deactivated => reload
        engine_now = self.engine_class.objects.get(
            id=self.engine_id,
            is_active=True)
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

            self._engine = self.engine_class.objects.load(
                id=self.engine_id,
                is_active=True)
        return self._engine

    def run(self, *args, **kwargs):
        return self.engine.tasks(*args, **kwargs)


class PaperEngineTask(EngineTask):
    engine_class = PaperEngine


class ThreadEngineTask(EngineTask):
    engine_class = ThreadEngine
