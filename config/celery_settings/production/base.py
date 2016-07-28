# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery_settings.common import *
from config.utils import get_dns_name_based_on_role


EMAIL_BACKEND = 'django_mailgun.MailgunBackend'
MAILGUN_ACCESS_KEY = env.str('MAILGUN_KEY')
MAILGUN_SERVER_NAME = 'mg.etalia.io'

# Celery
BROKER_URL = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=get_dns_name_based_on_role('master'),
)
CELERY_RESULT_BACKEND = 'amqp://{username}:{password}@{host}:5672//'.format(
    username=env.str('RABBITMQ_USERNAME'),
    password=env.str('RABBITMQ_PASSWORD'),
    host=get_dns_name_based_on_role('master'),
)