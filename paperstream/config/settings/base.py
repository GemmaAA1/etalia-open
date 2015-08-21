"""
Django settings for paperstream project.

Generated by 'django-admin startproject' using Django 1.8.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.8/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
from unipath import Path
from kombu import Queue, Exchange
from celery.schedules import crontab

ROOT_DIR = Path(__file__).ancestor(4)  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.child('paperstream')

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.8/howto/deployment/checklist/

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True
ALLOWED_HOSTS = []


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
    'core',
    'library',
    'populate',
    'consumers',
    'users',
    'feeds',
    'nlp',
    # 'comments',
    # 'networks',
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
    # 'default': {
    #     'ENGINE': 'django.db.backends.sqlite3',
    #     'NAME': os.path.join(ROOT_DIR, '../database/db.sqlite3'),
    # }
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
    str(APPS_DIR.child('templates')),
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
# STATIC_ROOT = APPS_DIR.child('static')

# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-url
STATIC_URL = '/static/'

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#std:setting-STATICFILES_DIRS
STATICFILES_DIRS = (
    str(APPS_DIR.child('static')),
)

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#staticfiles-finders
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)


# MEDIA CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-root
MEDIA_ROOT = APPS_DIR.child('media')

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
    'users.backends.mendeley.CustomMendeleyOAuth2',
    'users.backends.zotero.CustomZoteroOAuth',
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
SOCIAL_AUTH_EMAIL_VALIDATION_FUNCTION = 'users.mail.send_validation'
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
    'users.pipeline.require_primary',
    'social.pipeline.mail.mail_validation',
    'social.pipeline.user.create_user',
    'social.pipeline.social_auth.associate_user',
    'social.pipeline.social_auth.load_extra_data',
    'social.pipeline.user.user_details',
    'social.pipeline.debug.debug',
    # 'users.pipeline.update_user_lib',
    'users.pipeline.init_user',
    'users.pipeline.require_affiliation',
)

# TODO:
# ****
# API KEY TO MOVE IN ENV VARIABLE LATER
# ****

# Mendeley
SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_KEY = '1678'
SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_SECRET = 'caOrLU0DqOUC4wdD'
# Zotero
SOCIAL_AUTH_CUSTOM_ZOTERO_KEY = 'a7ecbff3d0bbe59abc4b'
SOCIAL_AUTH_CUSTOM_ZOTERO_SECRET = 'c5d0c178d9196e62bdbf'


# DISQUS
DISQUS_API_KEY = ''
DISQUS_WEBSITE_SHORTNAME = 'paperstream'
# TODO: Move to var env
DISQUS_SECRET_KEY = 'eMWsm6qeNkDHzdvLViScWPldyDVnmvAz4U79YjsCelOu58XnRPelrUimqTrhGrRw'
DISQUS_PUBLIC_KEY = 'w2W0iBEJwGE49PjupwQxDnfzC9ayliEvctiGwbmVb63uHIXNTZLgreJDNRvvBOap'

# CONSUMER CONFIGURATION
# ------------------------------------------------------------------------------

# Minimum delay between same journal is consumed twice (in period of scheduler.
# ie if period is a day, then unit is in days)
CONS_MIN_DELAY = 0
CONS_MAX_DELAY = 7

# In days, how many day in the past to look at when initializing database
CONS_INIT_PAST = 365


# NLP APP
# ------------------------------------------------------------------------------
NLP_CHUNK_SIZE = 10000
NLP_DATA_PATH = os.path.join(str(APPS_DIR), 'nlp', 'data')
NLP_DOC2VEC_PATH = os.path.join(str(APPS_DIR), 'nlp', 'mods')
NLP_LSH_PATH = os.path.join(str(APPS_DIR), 'nlp', 'lshfs')
NLP_MAX_VECTOR_SIZE = 300
NLP_MAX_KNN_NEIGHBORS = 10

FEED_JOURNAL_VECTOR_RATIO = 0.2

# FEED APP
# ------------------------------------------------------------------------------
FEEDS_SCORE_KEEP_N_PAPERS = 100
FEEDS_DISPLAY_N_PAPERS = 50


# CELERY
# ------------------------------------------------------------------------------

BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'
# BROKER_URL = 'pyamqp://'
# CELERY_RESULT_BACKEND = 'pyamqp://'
CELERY_ACCEPT_CONTENT = ['json', 'pickle']
CELERY_TASK_RESULT_EXPIRES = 60  # in seconds

CELERY_DEFAULT_QUEUE = 'default'
# embed_exchange = Exchange('embed', type='topic')
# consumer_exchange = Exchange('consumer', type='topic')
CELERY_QUEUES = (
    Queue('default', routing_key='task.#'),
    Queue('dbow', routing_key='dbow.#'),
    Queue('consumers', routing_key='consumers.#'),
)
CELERY_DEFAULT_EXCHANGE = 'tasks'
CELERY_DEFAULT_EXCHANGE_TYPE = 'topic'
CELERY_DEFAULT_ROUTING_KEY = 'task.default'

CELERY_ROUTES = {
    'nlp.tasks.dbow_embed_paper': {
        'queue': 'dbow',
        'routing_key': 'dbow.embed',
    },
    'nlp.tasks.dbow_lsh': {
        'queue': 'dbow',
        'routing_key': 'dbow.lsh',
    },
    'consumers.tasks.pubmed_run_all': {
        'queue': 'consumers',
        'routing_key': 'consumers.pubmed',
    },
    'consumers.tasks.arxiv_run_all': {
        'queue': 'consumers',
        'routing_key': 'consumers.arxiv',
    },
    'consumers.tasks.elsevier_run_all': {
        'queue': 'consumers',
        'routing_key': 'consumers.elsevier',
    },
}

CELERYBEAT_SCHEDULE = {
    'pubmed-once-a-day': {
        'task': 'consumers.tasks.pubmed_run_all',
        'schedule': crontab(minute=0, hour=0),  # daily at midnight
    },
    'arxiv-once-a-day': {
        'task': 'consumers.tasks.arxiv_run_all',
        'schedule': crontab(minute=0, hour=12),  # daily at 12pm
    },
    'elsevier-once-a-day': {
        'task': 'consumers.tasks.elsevier_run_all',
        'schedule': crontab(minute=0, hour=18),  # daily at 6pm
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
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'paperstream.log'),
            'formatter': 'verbose'
        },
        'populate': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'populate.log'),
            'formatter': 'verbose'
        },
        'nlp': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'nlp.log'),
            'formatter': 'verbose'
        },
        'feeds': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'feeds.log'),
            'formatter': 'verbose'
        },
        'users': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'users.log'),
            'formatter': 'verbose'
        },
        'consumers': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'consumers.log'),
            'formatter': 'verbose'
        },
        'celery': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(ROOT_DIR.child('logs'), 'celery.log'),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 50,  # 100 mb
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
        #     'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
        # },
        'populate': {
            'handlers': ['console', 'populate'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'consumers': {
            'handlers': ['console', 'consumers'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'nlp': {
            'handlers': ['console', 'nlp'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'feeds': {
            'handlers': ['console', 'feeds'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
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