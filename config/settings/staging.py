# -*- coding: utf-8 -*-
from .common import *

CONFIG_FILE = __file__

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['alpha-u6VcayvtcI.pubstream.io']

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

# S3 NLP buckets
NLP_DATA_BUCKET_NAME = 'pubstream-staging-nlp-data'
NLP_MODELS_BUCKET_NAME = 'pubstream-staging-nlp-models'
NLP_MS_BUCKET_NAME = 'pubstream-staging-nlp-ms'

# Static asset configuration
ROOT_DIR_m1 = ROOT_DIR.path() - 1
STATIC_ROOT = str(ROOT_DIR_m1.path('static'))
STATIC_URL = '/static/'

STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

# EMAIL backend
EMAIL_BACKEND = 'django_mailgun.MailgunBackend'
MAILGUN_ACCESS_KEY = env.str('MAILGUN_KEY')
MAILGUN_SERVER_NAME = 'mg.etalia.io'

# Invite mode switch
INVITE_MODE = False