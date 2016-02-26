# -*- coding: utf-8 -*-
from .common import *

CONFIG_FILE = __file__

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['www.etalia.io', 'etalia.io', 'alpha.etalia.io']

# Admins received email of the full exception information when DEGUB=False
ADMINS = [('Nicolas', 'nicolas.pannetier@gmail.com'),
          ]

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': env.str('AWS_RDS_DB_NAME'),
        'USER': env.str('AWS_RDS_USERNAME'),
        'PASSWORD': env.str('AWS_RDS_PASSWORD'),
        'HOST': env.str('AWS_RDS_HOSTNAME'),
        'PORT': env.str('AWS_RDS_PORT'),
    }
}

# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.sqlite3',
#         'NAME': str(ROOT_DIR.path('db.sqlite3')),
#     }
# }

# S3 NLP buckets
NLP_DATA_BUCKET_NAME = 'etalia-production-nlp-data'
NLP_MODELS_BUCKET_NAME = 'etalia-production-nlp-models'
NLP_MS_BUCKET_NAME = 'etalia-production-nlp-ms'

# Static asset configuration
ROOT_DIR_m1 = ROOT_DIR.path() - 1
STATIC_ROOT = str(ROOT_DIR_m1.path('static'))
STATIC_URL = '/static/'


STATICFILES_DIRS = (
    str(APPS_DIR.path('static/compiled')),
)

# EMAIL backend
EMAIL_BACKEND = 'django_mailgun.MailgunBackend'
MAILGUN_ACCESS_KEY = env.str('MAILGUN_KEY')
MAILGUN_SERVER_NAME = 'mg.etalia.io'

INVITE_MODE = False
if INVITE_MODE:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': str((ROOT_DIR - 1).path('db').path('invite.sqlite3')),
        }
    }

