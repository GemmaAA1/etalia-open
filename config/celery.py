# -*- coding: utf-8 -*-
"""
Celery auto detect task script. Task defined with decorator celery_app
located in a app/task.py files are auto detected and available to workers

"""
from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
from paperstream.nlp.models import Model
from paperstream.nlp.tasks_class import EmbedPaperTask, MostSimilarTask
from django.conf import settings
from celery import bootsteps

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings.development')

celery_app = Celery('paperstream')

# Using a string here means the worker will not have to
# pickle the object when using Windows.
celery_app.config_from_object('django.conf:settings')
celery_app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)


# Registering Model and MostSimilar task
# -----------------------------------------------------------------------------
# Add user options to control for nlp / mostsimilar registration of tasks
from celery.bin import Option
celery_app.user_options['worker'].add(
    Option('--init', dest='init', default=None, help='init MostSimilar (ms) or Model based tasks (nlp)')
)

class NLPBootstep(bootsteps.Step):
    """Bootstep to register task with or without initializing data upload depending
    on 'init' argument"""

    def __init__(self, worker, init, **options):
        # register specific tasks
        if init:
            if 'nlp' in init:
                self.register_model_tasks(init=True)
            if 'ms' in init:
                self.register_mostsimilar_tasks(init=True)
        else:
            self.register_model_tasks()
            self.register_mostsimilar_tasks()

    @staticmethod
    def register_model_tasks(init=False):
        """Register Model tasks
        """
        model_names = Model.objects.all().values_list('name', flat=True)
        for model_name in model_names:
            cls = EmbedPaperTask(model_name=model_name, init=init)
            celery_app.task(cls, name='paperstream.nlp.tasks.{model_name}'.format(
                model_name=model_name))

    @staticmethod
    def register_mostsimilar_tasks(init=False):
        """Register MostSimilar tasks
        """
        model_names = Model.objects.all().values_list('name', flat=True)
        for model_name in model_names:
            cls = MostSimilarTask(model_name=model_name, init=init)
            celery_app.task(cls, name='paperstream.nlp.tasks.mostsimilar_{model_name}'.format(
                model_name=model_name))


celery_app.steps['worker'].add(NLPBootstep)


@celery_app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))
