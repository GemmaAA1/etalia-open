# -*- coding: utf-8 -*-

from .common import *

CONFIG_FILE = __file__

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = env.bool('DJANGO_DEBUG', default=True)
TEMPLATE_DEBUG = DEBUG

# DEBUG TOOLBAR
# ------------------------------------------------------------------------------
DEBUG_TOOLBAR_PATCH_SETTINGS = False
# INSTALLED_APPS += ('debug_toolbar', )
# MIDDLEWARE_CLASSES = ('debug_toolbar.middleware.DebugToolbarMiddleware', ) + \
#                      MIDDLEWARE_CLASSES
INTERNAL_IPS = ['127.0.0.1']

# APPS
INSTALLED_APPS += ('django_extensions', )


GRAPH_MODELS = {
  'all_applications': True,
  'group_models': True,
}


# CONSUMER
# ------------------------------------------------------------------------------
# In days, how many day in the past to look at when initializing database
CONS_INIT_PAST = 30

# Mail settings
# ------------------------------------------------------------------------------
EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_BACKEND = env('DJANGO_EMAIL_BACKEND',
                    default='django.core.mail.backends.console.EmailBackend')

# Static asset configuration
STATIC_ROOT = str(ROOT_DIR.path('staticfiles'))
STATIC_URL = '/static/'
# STATIC_URL = str(ROOT_DIR.path('staticfiles'))

STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

# AWS S3 Buckets
NLP_DATA_BUCKET_NAME = 'etalia-development-nlp-data'
NLP_MODELS_BUCKET_NAME = 'etalia-development-nlp-models'
NLP_MS_BUCKET_NAME = 'etalia-development-nlp-ms'
AWS_STORAGE_BUCKET_NAME = env('DJANGO_AWS_STORAGE_BUCKET_NAME')


INVITE_MODE = False
if INVITE_MODE:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': str((ROOT_DIR - 1).path('db').path('invite.sqlite3')),
        }
    }

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        'TIMEOUT': 60 * 60 * 24,    # 24 h
    },
    'files': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        'TIMEOUT': 60 * 60,         # 1 h
    }
}

#!!!!! WARNING, USE WITH CARE !!!!!!####
#                                      #
#      SWITCH TO PRODUCTION DATABASE   #
#                                      #
########################################

#
# NLP_DATA_BUCKET_NAME = 'etalia-production-nlp-data'
# NLP_MODELS_BUCKET_NAME = 'etalia-production-nlp-models'
# NLP_MS_BUCKET_NAME = 'etalia-production-nlp-ms'
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
