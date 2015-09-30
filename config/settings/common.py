# -*- coding: utf-8 -*-

"""
Django settings for paperstream project.

Generated by 'django-admin startproject' using Django 1.8.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.8/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
from __future__ import absolute_import, unicode_literals
import environ
from kombu import Queue, Exchange
from celery.schedules import crontab

ROOT_DIR = environ.Path(__file__) - 3  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.path('paperstream')

env = environ.Env()

SITE_ID = 1

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.8/howto/deployment/checklist/

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = env('DJANGO_DEBUG', default=False)
SECRET_KEY = env('DJANGO_SECRET_KEY', default='CHANGEME!!!')
ALLOWED_HOSTS = ['localhost', '127.0.0.1', 'localhost']


# APP CONFIGURATION
# ------------------------------------------------------------------------------
DJANGO_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.sites',
)


THIRD_PARTY_APPS = (
    'social.apps.django_app.default',
    'disqus',
)

LOCAL_APPS = (
    'paperstream.core',
    'paperstream.library',
    'paperstream.populate',
    'paperstream.consumers',
    'paperstream.nlp',
    'paperstream.users',
    'paperstream.feeds',
    'paperstream.altmetric',
    # 'functional_tests',
)

INSTALLED_APPS = DJANGO_APPS + THIRD_PARTY_APPS + LOCAL_APPS

# MIDDLEWARE CONFIGURATION
# ------------------------------------------------------------------------------
MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django.middleware.security.SecurityMiddleware',
)


# DATABASE CONFIGURATION
# ------------------------------------------------------------------------------
# https://docs.djangoproject.com/en/1.8/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'paperstream',
        'USER': 'nicolaspannetier',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': '',
        'ATOMIC_REQUESTS': False,
        # NB: True conflicts with the use of python-social-auth (whose entire
        # pipeline is atomic while celery needs to know user during the pipeline
        # authentication process TODO: find a fix ?
    }
}

# TEMPLATE CONFIGURATION
# ------------------------------------------------------------------------------

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.contrib.auth.context_processors.auth',
    'django.core.context_processors.debug',
    'django.core.context_processors.i18n',
    'django.core.context_processors.media',
    'django.core.context_processors.static',
    'django.core.context_processors.tz',
    'django.contrib.messages.context_processors.messages',
    'django.core.context_processors.request',
    # Your stuff: custom template context processors go here
    'social.apps.django_app.context_processors.backends',
)

# See: https://docs.djangoproject.com/en/dev/ref/settings/#template-dirs
TEMPLATE_DIRS = (
    str(APPS_DIR.path('templates')),
)

TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
)

ITEMS_PER_PAGE = 15
NUMBER_OF_NEIGHBORS = 5

# GENERAL CONFIGURATION
# ------------------------------------------------------------------------------
# https://docs.djangoproject.com/en/1.8/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# STATIC FILE CONFIGURATION
#  https://docs.djangoproject.com/en/1.8/howto/static-files/
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-root
# use to serve static file in production by collecting static files in root
STATIC_ROOT = str(ROOT_DIR.path('staticfiles'))

# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-url
STATIC_URL = '/static/'

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#std:setting-STATICFILES_DIRS
STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#staticfiles-finders
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)


# MEDIA CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-root
MEDIA_ROOT = str(APPS_DIR.path('media'))

# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-url
MEDIA_URL = '/media/'


# URL Configuration
# ------------------------------------------------------------------------------
ROOT_URLCONF = 'config.urls'

# See: https://docs.djangoproject.com/en/dev/ref/settings/#wsgi-application
WSGI_APPLICATION = 'config.wsgi.application'


# AUTHENTICATION CONFIGURATION
# ------------------------------------------------------------------------------
AUTHENTICATION_BACKENDS = (
    'paperstream.users.backends.mendeley.CustomMendeleyOAuth2',
    'paperstream.users.backends.zotero.CustomZoteroOAuth',
    'social.backends.email.EmailAuth',
    'django.contrib.auth.backends.ModelBackend',
    # 'allauth.account.auth_backends.AuthenticationBackend',
)

AUTH_USER_MODEL = 'users.User'

# LOGIN_URL = '/login/'
LOGIN_URL = '/'
LOGIN_REDIRECT_URL = '/'
URL_PATH = ''
SOCIAL_AUTH_STRATEGY = 'social.strategies.django_strategy.DjangoStrategy'
SOCIAL_AUTH_STORAGE = 'social.apps.django_app.default.models.DjangoStorage'
SOCIAL_AUTH_GOOGLE_OAUTH_SCOPE = [
    'https://www.googleapis.com/auth/drive',
    'https://www.googleapis.com/auth/userinfo.profile'
]
# SOCIAL_AUTH_EMAIL_FORM_URL = '/signup-email'
SOCIAL_AUTH_EMAIL_FORM_HTML = 'email_signup.html'
SOCIAL_AUTH_EMAIL_VALIDATION_FUNCTION = 'paperstream.users.mail.send_validation'
SOCIAL_AUTH_EMAIL_VALIDATION_URL = 'user:validation-sent'
# SOCIAL_AUTH_USERNAME_FORM_URL = '/signup-username'
SOCIAL_AUTH_USERNAME_FORM_HTML = 'username_signup.html'
SOCIAL_AUTH_USERNAME_IS_FULL_EMAIL = True

EMAIL_FROM = 'noreply@paperstream.io'

SOCIAL_AUTH_PIPELINE = (
    'social.pipeline.social_auth.social_details',
    'social.pipeline.social_auth.social_uid',
    'social.pipeline.social_auth.auth_allowed',
    'social.pipeline.social_auth.social_user',
    'social.pipeline.user.get_username',
    'paperstream.users.pipeline.require_primary',
    # 'social.pipeline.mail.mail_validation',
    'social.pipeline.user.create_user',
    'social.pipeline.social_auth.associate_user',
    'social.pipeline.social_auth.load_extra_data',
    'social.pipeline.user.user_details',
    'social.pipeline.debug.debug',
    'paperstream.users.pipeline.init_user',
    # 'paperstream.users.pipeline.require_affiliation',
)

# Mendeley Keys
SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_KEY = \
    env('SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_KEY')
SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_SECRET = \
    env('SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_SECRET')

# Zotero Keys
SOCIAL_AUTH_CUSTOM_ZOTERO_KEY = env('SOCIAL_AUTH_CUSTOM_ZOTERO_KEY')
SOCIAL_AUTH_CUSTOM_ZOTERO_SECRET = env('SOCIAL_AUTH_CUSTOM_ZOTERO_SECRET')


# DISQUS
# ------------------------------------------------------------------------------
DISQUS_WEBSITE_SHORTNAME = env('DISQUS_WEBSITE_SHORTNAME')
DISQUS_PUBLIC_KEY = env('DISQUS_PUBLIC_KEY')
DISQUS_SECRET_KEY = env('DISQUS_SECRET_KEY')


# AWS S3
# ------------------------------------------------------------------------------
DJANGO_AWS_ACCESS_KEY_ID = env('DJANGO_AWS_ACCESS_KEY_ID')
DJANGO_AWS_SECRET_ACCESS_KEY = env('DJANGO_AWS_SECRET_ACCESS_KEY')


# CONSUMER CONFIGURATION
# ------------------------------------------------------------------------------

# Minimum delay between same journal is consumed twice (in period of scheduler.
# ie if period is a day, then unit is in days)
CONS_MIN_DELAY = 0
CONS_MAX_DELAY = 7

# In days, how many day in the past to look at when initializing database
CONS_INIT_PAST = 365

CONSUMER_PUBMED_EMAIL = env('CONSUMER_PUBMED_EMAIL')
CONSUMER_ELSEVIER_API_KEY = env('CONSUMER_ELSEVIER_API_KEY')



# NLP APP
# ------------------------------------------------------------------------------
NLP_CHUNK_SIZE = 10000
NLP_DATA_PATH = str(ROOT_DIR.path('nlp_data', 'data'))
NLP_MODELS_PATH = str(ROOT_DIR.path('nlp_data', 'mods'))
NLP_MS_PATH = str(ROOT_DIR.path('nlp_data', 'ms'))
NLP_NLTK_DATA_PATH = str(APPS_DIR.path('nlp', 'nltk_data'))
NLP_MAX_VECTOR_SIZE = 300
NLP_MAX_KNN_NEIGHBORS = 10

NLP_TIME_LAPSE_CHOICES = (
    (7, '1 Week'),
    (30, '1 Month'),
    (60, '2 Months'),
    (365, '1 Year'),
    (-1, 'All'),
)

# Time in days for recomputing neighbors is accessed
NLP_NEIGHBORS_REFRESH_TIME_LAPSE = 7

# FEED APP
# ------------------------------------------------------------------------------
FEEDS_SCORE_KEEP_N_PAPERS = 100
FEEDS_DISPLAY_N_PAPERS = 50
FEED_JOURNAL_VECTOR_RATIO = 0.2
# Number of neighbors from seed paper
FEED_K_NEIGHBORS = 10

# ALTMETRIC APP
# ------------------------------------------------------------------------------
ALTMETRIC_API_KEY = env('ALTMETRIC_API_KEY')
# slightly less than each second in a a day
ALTMETRIC_MAX_PAPERS_PER_PERIOD = 20 * 3600


# LANDING
# ------------------------------------------------------------------------------
LANDING_ACTIVE_PAPERS_NUMBER = 2
LANDING_ACTIVE_PAPERS_TIME_IN_DAYS = 30

# CELERY
# ------------------------------------------------------------------------------

BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'
CELERY_ACCEPT_CONTENT = ['json', 'pickle']
CELERY_TASK_RESULT_EXPIRES = 60  # in seconds

CELERY_DEFAULT_QUEUE = 'default'
# embed_exchange = Exchange('embed', type='topic')
# consumer_exchange = Exchange('consumer', type='topic')
CELERY_QUEUES = (
    Queue('default', routing_key='default.#'),
    Queue('nlp', routing_key='nlp.#'),
    Queue('mostsimilar', routing_key='mostsimilar.#'),
    Queue('consumers', routing_key='consumers.#'),
    Queue('altmetric', routing_key='altmetric.#'),
)
CELERY_DEFAULT_EXCHANGE = 'tasks'
CELERY_DEFAULT_EXCHANGE_TYPE = 'topic'
CELERY_DEFAULT_ROUTING_KEY = 'default'

CELERY_ROUTES = ('config.routers.MyRouter', )

CELERYBEAT_SCHEDULE = {
    'update-altmetric': {
        'task': 'paperstream.altmetric.tasks.update_altmetric_periodic',
        'schedule': crontab(minute=0, hour=0, day_of_week='*/2'),  # every 2 days at UTC+0
    },
    'update-ms-all': {
        'task': 'paperstream.nlp.tasks.mostsimilar_update_all',
        'schedule': crontab(minute=0, hour=0),  # every day at UTC+0
    },
    'pubmed-once-a-day': {
        'task': 'paperstream.consumers.tasks.pubmed_run_all',
        'schedule': crontab(minute=0, hour=6),  # daily at UCT+6
    },
    'arxiv-once-a-day': {
        'task': 'paperstream.consumers.tasks.arxiv_run_all',
        'schedule': crontab(minute=0, hour=12),  # daily at UTC+12
    },
    'elsevier-once-a-day': {
        'task': 'paperstream.consumers.tasks.elsevier_run_all',
        'schedule': crontab(minute=0, hour=18),  # daily at UTC+18
    },
}

# LOGGING CONFIGURATION
# ------------------------------------------------------------------------------
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': "[%(asctime)s.%(msecs)03d] %(levelname)s (%(name)s:%(funcName)s) %(message)s",
            'datefmt': "%Y-%m-%d %H:%M:%S",
        },
        'simple': {
            'format': '%(levelname)s %(module)s %(message)s'
        },
    },
    'handlers': {
        'null': {
            'level': 'DEBUG',
            'class': 'logging.NullHandler',
        },
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },
        'file': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'paperstream.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'populate': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'populate.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'nlp': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'nlp.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'feeds': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'feeds.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'users': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'users.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'consumers': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'consumers.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'celery': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(ROOT_DIR.path('logs', 'celery.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
    },
    'loggers': {
        # 'django': {
        #     'handlers': ['null'],
        #     'propagate': True,
        #     'level': 'INFO',
        # },
        # 'paperstream': {
        #     'handlers': ['console', 'file'],
        #     'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
        # },
        'populate': {
            'handlers': ['console', 'populate'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'consumers': {
            'handlers': ['console', 'consumers'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'nlp': {
            'handlers': ['console', 'nlp'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'feeds': {
            'handlers': ['console', 'feeds'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'users': {
            'handlers': ['console', 'users'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'celery': {
            'handlers': ['celery', 'console'],
            'level': 'DEBUG',
            'propagate': False,
        },
    }
}