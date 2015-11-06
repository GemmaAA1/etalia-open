# -*- coding: utf-8 -*-
from .common import *

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['www.pubstream.io']

# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.postgresql_psycopg2',
#         'NAME': env.str('AWS_RDS_DB_NAME'),
#         'USER': env.str('AWS_RDS_USERNAME'),
#         'PASSWORD': env.str('AWS_RDS_PASSWORD'),
#         'HOST': env.str('AWS_RDS_HOSTNAME'),
#         'PORT': env.str('AWS_RDS_PORT'),
#     }
# }

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': str(ROOT_DIR.path('db.sqlite3')),
    }
}

NLP_DATA_BUCKET_NAME = 'pubstream-production-nlp-data'
NLP_MODELS_BUCKET_NAME = 'pubstream-production-nlp-models'
NLP_MS_BUCKET_NAME = 'pubstream-production-nlp-ms'

INVITE_MODE = True
if INVITE_MODE:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': str((ROOT_DIR - 1).path('db').path('invite.sqlite3')),
        }
    }

