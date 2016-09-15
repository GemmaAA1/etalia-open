# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery_settings.common import *
from config.utils import get_private_ip_based_on_role


ANYMAIL = {
    "MAILGUN_API_KEY": env.str('MAILGUN_KEY'),
    "MAILGUN_SENDER_DOMAIN": 'mg.etalia.io'
}
EMAIL_BACKEND = 'anymail.backends.mailgun.MailgunBackend'
DEFAULT_FROM_EMAIL = 'contact@etalia.io'

# Celery
BROKER_URL = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=env.str('RABBITMQ_ELASTIC_IP'),
)
CELERY_RESULT_BACKEND = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=env.str('RABBITMQ_ELASTIC_IP'),
)