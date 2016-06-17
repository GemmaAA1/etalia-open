# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from .models import Model, PaperEngine, ThreadEngine
from .models import UserFingerprint
from .tasks_class import EmbedPaperTask, ThreadEngineTask, PaperEngineTask
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

# Register Model based tasks
#-------------------------------------------------------------------------------
try:
    model_name = Model.objects.get(is_active=True)
    # specific task name
    task_name ='etalia.nlp.tasks.nlp_dispatcher_{0}'.format(model_name)


    @app.task(base=EmbedPaperTask, bind=True, name=task_name, model_name=model_name)
    def nlp_dispatcher(self, *args, **kwargs):
        return self.model.tasks(*args, **kwargs)
except Model.DoesNotExist:
    pass

# Register PaperEngine based tasks
#-------------------------------------------------------------------------------
try:
    pe = PaperEngine.objects.get(is_active=True)
    # specific task name
    task_name ='etalia.nlp.tasks.pe_dispatcher_{0}'.format(pe.name)


    @app.task(base=PaperEngineTask, bind=True, name=task_name,
              engine_id=pe.id)
    def pe_dispatcher(self, *args, **kwargs):
        return self.engine.tasks(*args, **kwargs)
except PaperEngine.DoesNotExist:
    pass

# Register ThreadEngine based tasks
#-------------------------------------------------------------------------------
try:
    te = ThreadEngine.objects.get(is_active=True)
    # specific task name
    task_name ='etalia.nlp.tasks.te_dipatcher_{0}'.format(te.name)


    @app.task(base=ThreadEngineTask, bind=True, name=task_name,
              engine_id=te.id)
    def te_dispatcher(self, *args, **kwargs):
        return self.engine.tasks(*args, **kwargs)
except ThreadEngine.DoesNotExist:
    pass

# Update tasks
# ------------
def paperengine_update_all():
    pe_ids = PaperEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for peid in pe_ids:
        pe = PaperEngine.objects.load(id=peid, is_active=True)
        pe.update()


@app.task()
def paperengine_full_update_all():
    pe_ids = PaperEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for peid in pe_ids:
        pe = PaperEngine.objects.load(id=peid, is_active=True)
        pe.full_update()


@app.task()
def threadengine_update_all():
    te_ids = ThreadEngine.objects.filter(is_active=True).values_list('id', flat=True)
    for teid in te_ids:
        te = ThreadEngine.objects.load(id=teid, is_active=True)
        te.update()


@app.task()
def userfingerprints_update_all():
    ufps = UserFingerprint.objects.all()
    for ufp in ufps:
        ufp.update()


@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y
