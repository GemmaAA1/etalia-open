# -*- coding: utf-8 -*-
"""
Celery auto detect task script. Task defined with decorator celery_app
located in a app/task.py files are auto detected and available to workers

"""
from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
from paperstream.nlp.unreg_tasks import register_model_tasks, \
    register_mostsimilar_tasks
from django.conf import settings

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings.development')

celery_app = Celery('paperstream')

# Using a string here means the worker will not have to
# pickle the object when using Windows.
celery_app.config_from_object('django.conf:settings')
celery_app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)


# Add user options to control for nlp / mostsimilar registration of tasks
from celery import bootsteps
from celery.bin import Option
celery_app.user_options['worker'].add(
    Option('--role', dest='role', default=None, help='role for registering tasks')
)

class MyBootstep(bootsteps.Step):

    def __init__(self, worker, role=None, **options):

        # register role specific tasks
        if role:
            if 'nlp' in role:
                register_model_tasks()
            if 'ms' in role:
                register_mostsimilar_tasks()

celery_app.steps['worker'].add(MyBootstep)

@celery_app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))

