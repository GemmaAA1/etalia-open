# -*- coding: utf-8 -*-
from .common import *
from config.utils import get_dns_name_based_on_role

CONFIG_FILE = __file__

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['www.etalia.io', 'etalia.io', 'alpha.etalia.io',
                 'www.etalia.org', 'etalia.org']

SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

# Admins received email of the full exception information when DEGUB=False
ADMINS = [('Nicolas', 'webmaster@etalia.org'),
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
STATICFILES_STORAGE = 'django.contrib.staticfiles.storage.CachedStaticFilesStorage'

# EMAIL backend
ANYMAIL = {
    # "MAILGUN_API_KEY": env.str('MAILGUN_KEY'),
    # "MAILGUN_SENDER_DOMAIN": 'mg.etalia.org'
    "CUSTOMMAILGUN_API_KEY": env.str('MAILGUN_KEY'),
    "CUSTOMMAILGUN_SENDER_DOMAIN": 'mg.etalia.org'
}
# EMAIL_BACKEND = 'anymail.backends.mailgun.MailgunBackend'
EMAIL_BACKEND = 'etalia.core.emails.CustomMailgunBackend'
DEFAULT_FROM_EMAIL = 'contact@etalia.org'


# CACHE
# ------------------------------------------------------------------------------
SESSION_ENGINE = "django.contrib.sessions.backends.cached_db"

CACHE_FILE_DIR = str((ROOT_DIR-1).path('cache_files'))
# Create directory to store cache files if does not exists
if not os.path.exists(CACHE_FILE_DIR):
    os.mkdir(CACHE_FILE_DIR)


CACHES = {
    'default': {
        'BACKEND': 'django_redis.cache.RedisCache',
        'LOCATION': 'redis://{host}:6379/1'.format(host=env.str("REDIS_ELASTIC_IP")),
        'OPTIONS': {
            'CLIENT_CLASS': 'django_redis.client.DefaultClient',
            'IGNORE_EXCEPTIONS': True,
        },
        # 'TIMEOUT': 600,     # in seconds
    },
    'staticfiles': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'staticfiles-filehashes'
    },
    'files': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': CACHE_FILE_DIR,
        'TIMEOUT': 60 * 60,     # 1 h
    }
}


