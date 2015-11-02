# -*- coding: utf-8 -*-
from django.db import connection
from django.db.utils import ProgrammingError
from .celery import celery_app, register_model_tasks, register_mostsimilar_tasks

# To import tasks when running outside of celery
# Only register task if database has been created
try:
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM nlp_model')
    celery_app.loader.import_default_modules()
    register_model_tasks()
    register_mostsimilar_tasks()
except ProgrammingError:
    pass
