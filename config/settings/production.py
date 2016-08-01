# -*- coding: utf-8 -*-
from .common import *
from config.utils import get_dns_name_based_on_role

CONFIG_FILE = __file__

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['www.etalia.io', 'etalia.io', 'alpha.etalia.io']

SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

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
NLP_PE_BUCKET_NAME = 'etalia-production-nlp-pe'
NLP_TE_BUCKET_NAME = 'etalia-production-nlp-pe'

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


# CACHE
# ------------------------------------------------------------------------------
SESSION_ENGINE = "django.contrib.sessions.backends.cached_db"

CACHE_FILE_DIR = str((ROOT_DIR-1).path('cache_files'))
# Create directory to store cache files if does not exists
if not os.path.exists(CACHE_FILE_DIR):
    os.mkdir(CACHE_FILE_DIR)


CACHES = {
    'default': {
        'BACKEND': 'redis_cache.RedisCache',
        'LOCATION': '{host}:6379'.format(host=env.str("REDIS_ELASTIC_IP")),
        'TIMEOUT': 600,     # in seconds
    },
    'files': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': CACHE_FILE_DIR,
        'TIMEOUT': 60 * 60,     # 1 h
    }
}


