# -*- coding: utf-8 -*-

from .celery import celery_app

# To import tasks when running outside of celery
celery_app.loader.import_default_modules()


