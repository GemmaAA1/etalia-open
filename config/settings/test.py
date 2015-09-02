# -*- coding: utf-8 -*-

from .common import *

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = False
TEMPLATE_DEBUG = False

# Mail settings
# ------------------------------------------------------------------------------
EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_BACKEND = env.str('DJANGO_EMAIL_BACKEND',
                        'django.core.mail.backends.console.EmailBackend')

# CACHING
# ------------------------------------------------------------------------------
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': ''
    }
}

# CONSUMER CONFIGURATION
# ------------------------------------------------------------------------------

# NLP PATHS CHANGE
NLP_CHUNK_SIZE = 2
NLP_MAX_VECTOR_SIZE = 300
NLP_DATA_PATH = str(ROOT_DIR.path('nlp_data_test', 'data'))
NLP_DOC2VEC_PATH = str(ROOT_DIR.path('nlp_data_test', 'mods'))
NLP_LSH_PATH = str(ROOT_DIR.path('nlp_data_test', 'lshfs'))

# LOGGING
LOGGING = {}

# TRICKS
PASSWORD_HASHERS = (
    'django.contrib.auth.hashers.MD5PasswordHasher',
)

# CELERY
CELERY_ALWAYS_EAGER = True
CELERY_EAGER_PROPAGATES_EXCEPTIONS = True
BROKER_BACKEND = 'memory'
