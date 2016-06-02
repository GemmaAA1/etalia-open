# -*- coding: utf-8 -*-
"""
Celery auto detect task script. Task defined with decorator celery_app
located in a app/task.py files are auto detected and available to workers

"""
from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
from etalia.nlp.models import Model, PaperEngine, ThreadEngine
from etalia.nlp.tasks_class import EmbedPaperTask, PaperEngineTask, \
    ThreadEngineTask
from django.conf import settings
from celery import bootsteps
from celery.bin import Option

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings.development')

celery_app = Celery('etalia')

# If celery_app.conf  has a broker it has been configured at launch, otherwise
# load development settings
if not celery_app.conf['BROKER_URL']:
    if 'staging' in os.environ.get('DJANGO_SETTINGS_MODULE'):
        celery_app.config_from_object('config.celery_settings.staging.base')
    elif 'production' in os.environ.get('DJANGO_SETTINGS_MODULE'):
        celery_app.config_from_object('config.celery_settings.production.base')
    elif 'development' in os.environ.get('DJANGO_SETTINGS_MODULE'):
        celery_app.config_from_object('config.celery_settings.development')

celery_app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

# Registering Model and MostSimilar task and init them or not depending on user
# options --init
# -----------------------------------------------------------------------------
# Add user options to control for nlp / mostsimilar registration of tasks
celery_app.user_options['worker'].add(
    Option('--init', dest='init', default=None,
           help='init PaperEngine, ThreadEngine or Model based tasks')
)


class NLPBootstep(bootsteps.Step):
    """Bootstep to register task with or without initializing data upload depending
    on 'init' argument when launch from celery (cf supervisor conf file for examples)"""

    def __init__(self, worker, init, **options):
        # register specific tasks
        if init:
            init_nlp = True if 'nlp' in init else False
            init_pe = True if 'pe' in init else False
            init_te = True if 'te' in init else False
            register_model_tasks(init=init_nlp)
            register_paperengine_tasks(init=init_pe)
            register_threadengine_tasks(init=init_te)


def register_model_tasks(init=False):
    """Register Model tasks
    """
    model_names = Model.objects \
        .filter(is_active=True) \
        .values_list('name', flat=True)
    for model_name in model_names:
        cls = EmbedPaperTask(model_name=model_name, init=init)
        celery_app.task(cls, name='etalia.nlp.tasks.{model_name}'.format(
            model_name=model_name))


def register_paperengine_tasks(init=False):
    """Register PaperEngine tasks
    """
    pes = PaperEngine.objects \
        .filter(is_active=True)
    for pe in pes:
        cls = PaperEngineTask(engine_id=pe.id, init=init)
        celery_app.task(cls,
                        name='etalia.nlp.tasks.{name}'.format(name=pe.name))


def register_threadengine_tasks(init=False):
    """Register MostSimilar tasks
    """
    tes = ThreadEngine.objects \
        .filter(is_active=True)
    for te in tes:
        cls = ThreadEngineTask(engine_id=te.id, init=init)
        celery_app.task(cls,
                        name='etalia.nlp.tasks.{name}'.format(name=te.name)
                        )

celery_app.steps['worker'].add(NLPBootstep)


@celery_app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))
