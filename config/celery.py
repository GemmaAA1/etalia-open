# -*- coding: utf-8 -*-
"""
Celery auto detect task script. Task defined with decorator celery_app
located in a app/task.py files are auto detected and available to workers

"""
from __future__ import absolute_import, unicode_literals
import os
from celery import Celery, Task
from django.conf import settings
from celery import bootsteps
from celery.bin import Option

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings.development')

celery_app = Celery('etalia')

# If celery_app.conf  has a broker it has been configured at launch, otherwise
# load development settings
if not celery_app.conf.get('CELERY_BROKER_URL'):
    if 'production' in os.environ.get('DJANGO_SETTINGS_MODULE'):
        celery_app.config_from_object('config.celery_settings.production.base')
    elif 'development' in os.environ.get('DJANGO_SETTINGS_MODULE'):
        celery_app.config_from_object('config.celery_settings.development')

celery_app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

# -----------------------------------------------------------------------------
# # Add user options for load data at taks instantiation
# celery_app.user_options['worker'].add(
#     Option('--init', dest='init', default=None, help='Load data at task instantiation')
# )
#
#
# class InitArgs(bootsteps.Step):
#     """BootStep to warm up task dispatchers of type Model, PaperEngine or
#     ThreadEngine with data"""
#
#     def __init__(self, worker, init, **kwargs):
#         super(InitArgs, self).__init__(worker, **kwargs)
#         for k, task in worker.app.tasks.items():
#             if task.__name__.startswith('{0}'.format(init)):
#                 if init.startswith('nlp_dispatcher'):
#                     _ = task.model
#                 else:
#                     _ = task.engine
#
#
# celery_app.steps['worker'].add(InitArgs)

