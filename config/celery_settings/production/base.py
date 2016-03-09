# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery_settings.common import *

EMAIL_BACKEND = 'django_mailgun.MailgunBackend'
MAILGUN_ACCESS_KEY = env.str('MAILGUN_KEY')
MAILGUN_SERVER_NAME = 'mg.etalia.io'

# Celery
BROKER_URL = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=env.str('RABBITMQ_HOSTNAME'),
)
CELERY_RESULT_BACKEND = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=env.str('RABBITMQ_HOSTNAME'),
)
