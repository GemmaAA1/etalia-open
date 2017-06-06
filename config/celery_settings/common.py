# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
from kombu import Queue, Exchange
from kombu.common import Broadcast
from celery.schedules import crontab

ROOT_DIR = environ.Path(__file__) - 3  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.path('etalia')

env = environ.Env()

# EMAIl
CELERY_SEND_TASK_ERROR_EMAILS = True
SERVER_EMAIL = 'no-reply@etalia.org'
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
    Queue('beat', routing_key='beat.#'),
    Queue('test', routing_key='test.#'),
    Broadcast('broadcast_pe_tasks'),
    Broadcast('broadcast_te_tasks'),
)
CELERY_DEFAULT_EXCHANGE = 'tasks'
CELERY_DEFAULT_EXCHANGE_TYPE = 'topic'
CELERY_DEFAULT_ROUTING_KEY = 'default'

CELERY_ROUTES = ('config.routers.MyRouter', )

CELERYBEAT_SCHEDULE = {
    'update-altmetric': {
        'task': 'etalia.altmetric_app.tasks.update_altmetric_periodic',
        'schedule': crontab(minute=0, hour=0, day_of_week='mon,wed,fri'),
        'options': {'queue': 'beat'}
    },
    'update-userlib-all': {
        'task': 'etalia.users.tasks.userlib_update_all',
        'schedule': crontab(minute=0, hour=1),  # every day at UTC+1
        'options': {'queue': 'beat'}
    },
    'update-pe-all': {
        'task': 'etalia.nlp.tasks.paperengine_update_all',
        'schedule': crontab(minute=0, hour=2),  # every day at UTC+2
        'options': {'queue': 'beat'}
    },
    'update-te-all': {
        'task': 'etalia.nlp.tasks.threadengine_update_all',
        'schedule': crontab(minute=0, hour=3),  # every day at UTC+3
        'options': {'queue': 'beat'}
    },
    'update-all-userfingerprints': {
        'task': 'etalia.nlp.tasks.userfingerprints_update_all',
        'schedule': crontab(minute=0, hour=4, day_of_week='tue,thu,sat'),
        'options': {'queue': 'beat'}
    },
    'update-all-main-streams': {
        'task': 'etalia.feeds.tasks.update_all_main_streams',
        'schedule': crontab(minute=0, hour=5),  # every day at UTC+5
        'options': {'queue': 'beat'}
    },
    'update-all-main-trends': {
        'task': 'etalia.feeds.tasks.update_all_main_trends',
        'schedule': crontab(minute=0, hour=6),  # every day at UTC+6
        'options': {'queue': 'beat'}
    },
    'pubmed-once-a-day': {
        'task': 'etalia.consumers.tasks.pubmed_run_all',
        'schedule': crontab(minute=0, hour=7),  # daily at UTC+6
        'options': {'queue': 'beat'}
    },
    'arxiv-once-a-day': {
        'task': 'etalia.consumers.tasks.arxiv_run_all',
        'schedule': crontab(minute=0, hour=12),  # daily at UTC+12
        'options': {'queue': 'beat'}
    },
    'elsevier-once-a-day': {
        'task': 'etalia.consumers.tasks.elsevier_run_all',
        'schedule': crontab(minute=0, hour=18),  # daily at UTC+18
        'options': {'queue': 'beat'}
    },
    'springer-once-a-day': {
        'task': 'etalia.consumers.tasks.springer_run_all',
        'schedule': crontab(minute=0, hour=19),  # daily at UTC+19
        'options': {'queue': 'beat'}
    },
    'biorxiv-once-a-day': {
        'task': 'etalia.consumers.tasks.biorxiv_run_all',
        'schedule': crontab(minute=0, hour=20),  # daily at UTC+20
        'options': {'queue': 'beat'}
    },
    'populate-pubpeer': {
        'task': 'etalia.consumers.tasks.populate_pubpeer',
        'schedule': crontab(minute=0, hour=21),  # daily at UTC+21
        'options': {'queue': 'beat'}
    },
    'consolidate-library': {
        'task': 'etalia.consumers.tasks.consolidate_library',
        'schedule': crontab(minute=0, hour=20),  # daily at UTC-4
        'options': {'queue': 'beat'}
    },
    'emails-recommendations': {
        'task': 'etalia.users.tasks.send_recommendation_emails_on_wed_11am',
        'schedule': crontab(minute=0, hour=0, day_of_week='mon'),
        'options': {'queue': 'beat'}
    },
}