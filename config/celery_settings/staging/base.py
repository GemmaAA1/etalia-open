# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery_settings.common import *

# Celery
BROKER_URL = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=env.str('RABBITMQ_HOSTNAME'),
)
CELERY_RESULT_BACKEND = 'amqp://'