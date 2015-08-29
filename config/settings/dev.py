# -*- coding: utf-8 -*-

from .common import *

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = env.bool('DJANGO_DEBUG', default=True)
TEMPLATE_DEBUG = DEBUG


# Mail settings
# ------------------------------------------------------------------------------
EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_BACKEND = env('DJANGO_EMAIL_BACKEND',
                    default='django.core.mail.backends.console.EmailBackend')

# CACHING
# ------------------------------------------------------------------------------
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': ''
    }
}

# Static asset configuration
STATIC_ROOT = str(ROOT_DIR.path('staticfiles'))
STATIC_URL = '/static/'

STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

# STATICFILES_STORAGE = 'whitenoise.django.GzipManifestStaticFilesStorage'

# CELERY
# ------------------------------------------------------------------------------
# In development, all tasks will be executed locally by blocking until the task
# returns
CELERY_ALWAYS_EAGER = True

