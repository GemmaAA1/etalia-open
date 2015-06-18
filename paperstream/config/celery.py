"""
Celery auto detect task script. Task defined with decorator celery_app
located in a app/task.py files are auto detected and available to workers

launch with:
celery -A config worker --loglevel=info
"""

import os
from celery import Celery
from django.conf import settings

# set the default Django settings module for the 'celery' program.
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.local")

celery_app = Celery('paperstream')

# Using a string here means the worker will not have to
# pickle the object when using Windows.
celery_app.config_from_object('django.conf:settings')
celery_app.autodiscover_tasks(['library', 'consumers', 'users', 'feeds'])

@celery_app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))