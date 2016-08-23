# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
env = environ.Env()

ELASTIC_IP_MAPPING = {
    'redis': {
        'ip': env.str('REDIS_ELASTIC_IP'),
        'allocation_id': env.str('REDIS_ELASTIC_ALLOCATION_ID')
    },
    'master': {
        'ip': env.str('RABBITMQ_ELASTIC_IP'),
        'allocation_id': env.str('RABBITMQ_ELASTIC_ALLOCATION_ID')
    },
    'web': {
        'ip': env.str('WEB_ELASTIC_IP'),
        'allocation_id': env.str('WEB_ELASTIC_ALLOCATION_ID')
    }
}
