from .base import *
from core.utils import get_env_variable

# APPLICATION
INSTALLED_APPS += ('django_nose', )

# DEBUG
# ------------------------------------------------------------------------------
DEBUG = False
TEMPLATE_DEBUG = False


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


# LOGGING
LOGGING = {}

# TRICKS
PASSWORD_HASHERS = (
    'django.contrib.auth.hashers.MD5PasswordHasher',
)

# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.postgresql_psycopg2',
#         'NAME': 'paperstream',
#         'USER': 'nicolaspannetier',
#         'PASSWORD': '',
#         'HOST': 'localhost',
#         'PORT': '',
#         'ATOMIC_REQUESTS': False,
#         # NB: True conflicts with the use of python-social-auth (whose entire
#         # pipeline is atomic while celery needs to know user during the pipeline
#         # authentication process TODO: find a fix ?
#     }
# }


# CELERY
CELERY_ALWAYS_EAGER = True
CELERY_EAGER_PROPAGATES_EXCEPTIONS = True
BROKER_BACKEND = 'memory'

# Django Nose
TEST_RUNNER = 'django_nose.NoseTestSuiteRunner'