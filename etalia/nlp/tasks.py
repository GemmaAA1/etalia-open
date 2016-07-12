# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model
from django.db.models import Count, F
from .models import Model, PaperEngine, ThreadEngine
from .models import UserFingerprint
from .tasks_class import EmbedPaperTask, ThreadEngineTask, PaperEngineTask
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

User = get_user_model()

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
    task_name ='etalia.nlp.tasks.te_dispatcher_{0}'.format(te.name)


    @app.task(base=ThreadEngineTask, bind=True, name=task_name,
              engine_id=te.id)
    def te_dispatcher(self, *args, **kwargs):
        return self.engine.tasks(*args, **kwargs)
except ThreadEngine.DoesNotExist:
    pass


# Update tasks
# ------------
@app.task()
def paperengine_update_all():
    pe = PaperEngine.objects.load(is_active=True)
    pe.update()
    # pe_dispatcher.delay('update')


@app.task()
def paperengine_full_update_all():
    pe = PaperEngine.objects.load(is_active=True)
    pe.full_update()
    # pe_dispatcher.delay('full_update')


@app.task()
def threadengine_update_all():
    te = ThreadEngine.objects.load(is_active=True)
    te.update()
    # te_dispatcher.delay('update')


@app.task()
def userfingerprints_update_all():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_userfingerprint.delay(user_pk)


@app.task()
def update_userfingerprint(user_pk, name='main'):
    f, _ = UserFingerprint.objects.get_or_create(user_id=user_pk, name=name)
    f.update()
    return user_pk


@app.task()
def add_nlp(x, y):
    """dummy task"""
    logger.info("--> Processing task add")
    return x + y
