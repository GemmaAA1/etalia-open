# -*- coding: utf-8 -*-

from .celery import celery_app, register_model_tasks, register_mostsimilar_tasks

# To import tasks when running outside of celery
celery_app.loader.import_default_modules()
register_model_tasks()
register_mostsimilar_tasks()
