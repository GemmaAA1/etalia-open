# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
from kombu import Queue, Exchange
from celery.schedules import crontab

ROOT_DIR = environ.Path(__file__) - 3  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.path('etalia')

env = environ.Env()

# EMAIl
CELERY_SEND_TASK_ERROR_EMAILS = True
SERVER_EMAIL = 'no-reply@etalia.io'
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
ADMINS = [('Nicolas', 'nicolas.pannetier@gmail.com'),
          ]
# BROKER
BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'
CELERY_ACCEPT_CONTENT = ['json', 'pickle']
CELERY_TASK_RESULT_EXPIRES = 30  # in seconds

CELERY_DEFAULT_QUEUE = 'default'
CELERY_QUEUES = (
    Queue('default', routing_key='default.#'),
    Queue('nlp', routing_key='nlp.#'),
    Queue('pe', routing_key='pe.#'),
    Queue('te', routing_key='te.#'),
    Queue('update_engines', routing_key='update.#'),
    Queue('feed', routing_key='feed.#'),
    Queue('consumers', routing_key='consumers.#'),
    Queue('altmetric', routing_key='altmetric.#'),
    Queue('library', routing_key='library.#'),
    Queue('test', routing_key='test.#'),
)
CELERY_DEFAULT_EXCHANGE = 'tasks'
CELERY_DEFAULT_EXCHANGE_TYPE = 'topic'
CELERY_DEFAULT_ROUTING_KEY = 'default'

CELERY_ROUTES = ('config.routers.MyRouter', )

CELERYBEAT_SCHEDULE = {
    'update-altmetric': {
        'task': 'etalia.altmetric.tasks.update_altmetric_periodic',
        'schedule': crontab(minute=0, hour=0, day_of_week='*/2'),  # every 2 days at UTC+0
    },
    'update-userlib-all': {
        'task': 'etalia.users.tasks.userlib_update_all',
        'schedule': crontab(minute=0, hour=1),  # every day at UTC+1
    },
    'update-pe-all': {
        'task': 'etalia.nlp.tasks.paperengine_update_all',
        'schedule': crontab(minute=0, hour=2),  # every day at UTC+2
    },
    'update-te-all': {
        'task': 'etalia.nlp.tasks.threadengine_update_all',
        'schedule': crontab(minute=0, hour=3),  # every day at UTC+3
    },
    'update-all-userfingerprints': {
        'task': 'etalia.nlp.tasks.userfingerprints_update_all',
        'schedule': crontab(minute=0, hour=4, day_of_week='*/2'),  # every 2 days at UTC+4
    },
    'update-all-main-streams': {
        'task': 'etalia.feeds.tasks.update_all_main_streams',
        'schedule': crontab(minute=0, hour=5),  # every day at UTC+5
    },
    'update-all-main-trends': {
        'task': 'etalia.feeds.tasks.update_all_main_trends',
        'schedule': crontab(minute=0, hour=6),  # every day at UTC+6
    },
    'pubmed-once-a-day': {
        'task': 'etalia.consumers.tasks.pubmed_run_all',
        'schedule': crontab(minute=0, hour=7),  # daily at UTC+6
    },
    'arxiv-once-a-day': {
        'task': 'etalia.consumers.tasks.arxiv_run_all',
        'schedule': crontab(minute=0, hour=12),  # daily at UTC+12
    },
    'elsevier-once-a-day': {
        'task': 'etalia.consumers.tasks.elsevier_run_all',
        'schedule': crontab(minute=0, hour=18),  # daily at UTC+18
    },
    'consolidate-library': {
        'task': 'etalia.library.tasks.consolidate_library',
        'schedule': crontab(minute=0, hour=2),  # daily at UTC+2
    },
}