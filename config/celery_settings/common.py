# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
from kombu import Queue, Exchange
from celery.schedules import crontab

ROOT_DIR = environ.Path(__file__) - 3  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.path('paperstream')

env = environ.Env()

BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'
CELERY_ACCEPT_CONTENT = ['json', 'pickle']
CELERY_TASK_RESULT_EXPIRES = 5  # in seconds

CELERY_DEFAULT_QUEUE = 'default'
# embed_exchange = Exchange('embed', type='topic')
# consumer_exchange = Exchange('consumer', type='topic')
CELERY_QUEUES = (
    Queue('default', routing_key='default.#'),
    Queue('nlp', routing_key='nlp.#'),
    Queue('mostsimilar', routing_key='mostsimilar.#'),
    Queue('feed', routing_key='feed.#'),
    Queue('consumers', routing_key='consumers.#'),
    Queue('altmetric', routing_key='altmetric.#'),
)
CELERY_DEFAULT_EXCHANGE = 'tasks'
CELERY_DEFAULT_EXCHANGE_TYPE = 'topic'
CELERY_DEFAULT_ROUTING_KEY = 'default'

CELERY_ROUTES = ('config.routers.MyRouter', )

CELERYBEAT_SCHEDULE = {
    'update-library-stats': {
        'task': 'paperstream.library.tasks.update_stats',
        'schedule': crontab(minute=0, hour=0),  # every day at UTC+0
    },
    'update-altmetric': {
        'task': 'paperstream.altmetric.tasks.update_altmetric_periodic',
        'schedule': crontab(minute=0, hour=0, day_of_week='*/2'),  # every 2 days at UTC+0
    },
    'update-ms-all': {
        'task': 'paperstream.nlp.tasks.mostsimilar_update_all',
        'schedule': crontab(minute=0, hour=0),  # every day at UTC+0
    },
    'update-all-main-streams': {
        'task': 'paperstream.feeds.tasks.update_all_main_streams',
        'schedule': crontab(minute=0, hour=1),  # every day at UTC+1
    },
    'update-all-main-trends': {
        'task': 'paperstream.feeds.tasks.update_all_main_trends',
        'schedule': crontab(minute=0, hour=2),  # every day at UTC+2
    },
    'pubmed-once-a-day': {
        'task': 'paperstream.consumers.tasks.pubmed_run_all',
        'schedule': crontab(minute=0, hour=6),  # daily at UCT+6
    },
    'arxiv-once-a-day': {
        'task': 'paperstream.consumers.tasks.arxiv_run_all',
        'schedule': crontab(minute=0, hour=12),  # daily at UTC+12
    },
    'elsevier-once-a-day': {
        'task': 'paperstream.consumers.tasks.elsevier_run_all',
        'schedule': crontab(minute=0, hour=18),  # daily at UTC+18
    },
}