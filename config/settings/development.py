# -*- coding: utf-8 -*-

from .common import *

CONFIG_FILE = __file__
# CONFIG_FILE = 'production'

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = env.bool('DJANGO_DEBUG', default=True)
TEMPLATE_DEBUG = DEBUG

# DATABASE
# ------------------------------------------------------------------------------
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': os.environ.get('POSTGRES_DB', ''),
        'USER': os.environ.get('POSTGRES_USER', ''),
        'PASSWORD': os.environ.get('POSTGRES_PASSWORD', ''),
        'HOST': 'db',
        'PORT': 5432,
    },
}

# DEBUG TOOLBAR
# ------------------------------------------------------------------------------
DEBUG_TOOLBAR = False
if DEBUG_TOOLBAR:
    DEBUG_TOOLBAR_PATCH_SETTINGS = False
    INSTALLED_APPS += ('debug_toolbar', )
    MIDDLEWARE_CLASSES = ('debug_toolbar.middleware.DebugToolbarMiddleware', ) + \
                         MIDDLEWARE_CLASSES
    INTERNAL_IPS = ['127.0.0.1']

# APPS
# ------------------------------------------------------------------------------
# INSTALLED_APPS += ()


# CONSUMER
# ------------------------------------------------------------------------------
# In days, how many day in the past to look at when initializing database
CONSUMER_INIT_PAST = 7
CONSUMER_PUBPEER_INIT_PAST = 2

# Mail settings
# ------------------------------------------------------------------------------
DEFAULT_FROM_EMAIL = 'contact@etalia.io'

EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_BACKEND = env('DJANGO_EMAIL_BACKEND',
                    default='django.core.mail.backends.console.EmailBackend')
# For debug purposes only
ANYMAIL = {
    "CUSTOMMAILGUN_API_KEY": env.str('MAILGUN_KEY', default=''),
    "CUSTOMMAILGUN_SENDER_DOMAIN": 'mg.etalia.io'
}

# Static asset configuration
STATIC_ROOT = str(ROOT_DIR.path('staticfiles'))
STATIC_URL = '/static/'
STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

if 'production' in CONFIG_FILE:
    ROOT_DIR_m1 = ROOT_DIR.path() - 1
    STATIC_ROOT = str(ROOT_DIR_m1.path('static'))
    STATIC_URL = '/static/'


    STATICFILES_DIRS = (
        str(APPS_DIR.path('static/compiled')),
    )

# AWS S3 Buckets
NLP_DATA_BUCKET_NAME = 'etalia-development-nlp-data'
NLP_MODELS_BUCKET_NAME = 'etalia-development-nlp-models'
NLP_PE_BUCKET_NAME = 'etalia-development-nlp-pe'
NLP_TE_BUCKET_NAME = 'etalia-development-nlp-te'
AWS_STORAGE_BUCKET_NAME = env('AWS_STORAGE_BUCKET_NAME')

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        'TIMEOUT': 60 * 10,    # in seconds
    },
    'files': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        'TIMEOUT': 60 * 60,         # 1 h
    }
}


# Local redis cache for debug (redis-server must be on)
# CACHES = {
#     'default': {
#         'BACKEND': 'django_redis.cache.RedisCache',
#         'LOCATION': 'redis://{host}:6379/1'.format(host='127.0.0.1'),
#         'OPTIONS': {
#             'CLIENT_CLASS': 'django_redis.client.DefaultClient',
#             'IGNORE_EXCEPTIONS': True,
#         },
#         'TIMEOUT': 600,     # in seconds
#     },
# }
#     'default': {
#         'BACKEND': 'redis_cache.RedisCache',
#         'LOCATION': '{host}:6379'.format(host='localhost'),
#         'TIMEOUT': 300,     # 24 h
#     },
#     'files': {
#         'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
#         'TIMEOUT': 60 * 60,         # 1 h
#     }


#!!!!! WARNING, USE WITH CARE !!!!!!####
#                                      #
#      SWITCH TO PRODUCTION DATABASE   #
#                                      #
########################################


# NLP_DATA_BUCKET_NAME = 'etalia-production-nlp-data'
# NLP_MODELS_BUCKET_NAME = 'etalia-production-nlp-models'
# NLP_MS_BUCKET_NAME = 'etalia-production-nlp-pe'
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
#
# CACHE_FILE_DIR = str((ROOT_DIR-1).path('cache_files'))
# # Create directory to store cache files if does not exists
# if not os.path.exists(CACHE_FILE_DIR):
#     os.mkdir(CACHE_FILE_DIR)
#
# CACHES = {
#     'default': {
#         'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
#         'TIMEOUT': 60*5,     # 5 min
#     },
#     'files': {
#         'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
#         'LOCATION': CACHE_FILE_DIR,
#         'TIMEOUT': 60 * 60,     # 1 h
#     }
# }