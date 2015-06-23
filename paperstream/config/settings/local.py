from .base import *

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = get_env_variable('DJANGO_DEBUG', True)
TEMPLATE_DEBUG = DEBUG


# SECRET CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
# Note: This key only used for development and testing.
SECRET_KEY = get_env_variable('PAP_SECRET_KEY')


# Mail settings
# ------------------------------------------------------------------------------
EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_BACKEND = get_env_variable('DJANGO_EMAIL_BACKEND',
                    default='django.core.mail.backends.console.EmailBackend')

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

PUBMED_EMAIL = get_env_variable('PUBMED_EMAIL')
ELSEVIER_API_KEY = get_env_variable('ELSEVIER_API_KEY')