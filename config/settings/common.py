# -*- coding: utf-8 -*-

"""
Django settings for etalia project.

Generated by 'django-admin startproject' using Django 1.8.

For more information on this file, see
https://docs.djangoproject.com/en/1.8/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.8/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
from __future__ import absolute_import, unicode_literals
import environ
import os
from . import get_version

ROOT_DIR = environ.Path(__file__) - 3  # (/a/b/myfile.py - 3 = /)
APPS_DIR = ROOT_DIR.path('etalia')

env = environ.Env()

SITE_ID = 1

# Get app version from root __init__
VERSION = get_version(str(ROOT_DIR.path()))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.8/howto/deployment/checklist/

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = env.bool('DJANGO_DEBUG', default=False)
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
    'django.contrib.sitemaps',
)


THIRD_PARTY_APPS = (
    'social.apps.django_app.default',
    'avatar',
    'rest_framework',
    'rest_condition',
    'autofixture',
    'django_extensions',
    'rest_framework.authtoken',
    'crispy_forms',
    'easy_timezones',
    'anymail',
)

LOCAL_APPS = (
    'etalia.core',
    'etalia.users',
    'etalia.library',
    'etalia.threads',
    'etalia.nlp',
    'etalia.consumers',
    'etalia.feeds',
    'etalia.altmetric_app',
    'etalia.last_seen',
    'etalia.popovers',
    'etalia.usersession',
    'etalia.populate',
    'etalia.press',
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
    'easy_timezones.middleware.EasyTimezoneMiddleware',
    'etalia.last_seen.middleware.LastSeenMiddleware',
)


# DATABASE CONFIGURATION
# ------------------------------------------------------------------------------
# https://docs.djangoproject.com/en/1.8/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': '',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
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
    'social.apps.django_app.context_processors.backends',
    'django.core.context_processors.request',
    'etalia.core.context_processors.admin_context',
    'etalia.core.context_processors.user_update_check',
)

# See: https://docs.djangoproject.com/en/dev/ref/settings/#template-dirs
TEMPLATE_DIRS = (
    str(APPS_DIR.path('templates')),
    str(ROOT_DIR.path('campaigns/templates')),
)

TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
)


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


# Admin
# ------------------------------------------------------------------------------
# Admins received email of the full exception information when DEGUB=False
ADMINS = [('Nicolas', 'nicolas.pannetier@gmail.com'),
          ]

# EMAIL TEMPLATES
# ------------------------------------------------------------------------------
EMAIL_STATIC_BUCKET = 'https://s3-us-west-2.amazonaws.com/etalia-email-static/'
INVITE_EMAIL_TEMPLATE = 'emails/beta_invite_template-raw.html'
WELCOME_EMAIL_TEMPLATE = 'emails/welcome-raw.html'
INVITE_THREAD_EMAIL_TEMPLATE = 'emails/invite_thread-raw.html'


# RECOMMENDATION EMAILS
# ------------------------------------------------------------------------------

RECOMMENDATIONS_EMAILS_ON = True
PERIODIC_RECOMMENDATION_TEMPLATE = 'emails/periodic_recommendations-raw.html'
PERIODIC_RECOMMENDATION_NUMBER_PAPERS = 6


# AUTHENTICATION CONFIGURATION
# ------------------------------------------------------------------------------
AUTHENTICATION_BACKENDS = (
    'etalia.users.backends.mendeley.CustomMendeleyOAuth2',
    'etalia.users.backends.zotero.CustomZoteroOAuth',
    'etalia.users.backends.orcid.OrcidOAuth2',
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
SOCIAL_AUTH_EMAIL_VALIDATION_FUNCTION = 'etalia.users.mail.send_validation'
SOCIAL_AUTH_EMAIL_VALIDATION_URL = 'user:validation-sent'
# SOCIAL_AUTH_USERNAME_FORM_URL = '/signup-username'
SOCIAL_AUTH_USERNAME_FORM_HTML = 'username_signup.html'
SOCIAL_AUTH_USERNAME_IS_FULL_EMAIL = True

EMAIL_FROM = 'noreply@etalia.io'

SOCIAL_AUTH_PIPELINE = (
    'social.pipeline.social_auth.social_details',
    'social.pipeline.social_auth.social_uid',
    'social.pipeline.social_auth.auth_allowed',
    'social.pipeline.social_auth.social_user',
    'social.pipeline.user.get_username',
    'etalia.users.pipeline.require_primary',
    # 'social.pipeline.mail.mail_validation',
    'social.pipeline.user.create_user',
    'social.pipeline.social_auth.associate_user',
    'social.pipeline.social_auth.load_extra_data',
    'social.pipeline.user.user_details',
    # 'social.pipeline.debug.debug',
    'etalia.users.pipeline.create_details',
    'etalia.users.pipeline.update_usersession',
    'etalia.users.pipeline.send_email_of_new_signup',
    'etalia.users.pipeline.update_user',
    'etalia.users.pipeline.send_welcome_email_at_signup',
    # 'etalia.users.pipeline.require_affiliation',
)

# Mendeley Keys
SOCIAL_AUTH_MENDELEY_KEY = env('SOCIAL_AUTH_MENDELEY_KEY')
SOCIAL_AUTH_MENDELEY_SECRET = env('SOCIAL_AUTH_MENDELEY_SECRET')

# Zotero Keys
SOCIAL_AUTH_ZOTERO_KEY = env('SOCIAL_AUTH_ZOTERO_KEY')
SOCIAL_AUTH_ZOTERO_SECRET = env('SOCIAL_AUTH_ZOTERO_SECRET')
SOCIAL_AUTH_ZOTERO_AUTH_EXTRA_ARGUMENTS = {'write_access': '1'}

# ORCiD Keys
SOCIAL_AUTH_ORCID_KEY = 'APP-3S58TQE379Y00NXR'
SOCIAL_AUTH_ORCID_SECRET = '19826d6d-1648-4420-b101-11919a6cfce9'

# SESSION
# ------------------------------------------------------------------------------
SESSION_CACHE_ALIAS = 'default'


# LAST SEEN
# ------------------------------------------------------------------------------
LAST_SEEN_DEFAULT_MODULE = 'default'
LAST_SEEN_INTERVAL = 2 * 60 * 60


# AWS S3
# ------------------------------------------------------------------------------
AWS_ACCESS_KEY_ID = env('AWS_ACCESS_KEY_ID')
AWS_SECRET_ACCESS_KEY = env('AWS_SECRET_ACCESS_KEY')
AWS_STORAGE_BUCKET_NAME = env('AWS_STORAGE_BUCKET_NAME')
DEFAULT_FILE_STORAGE = 'storages.backends.s3boto.S3BotoStorage'

# CONSUMER CONFIGURATION
# ------------------------------------------------------------------------------

# Minimum delay between same journal is consumed twice (in period of scheduler.
# ie if period is a day, then unit is in days)
CONSUMER_MIN_DELAY = 0
CONSUMER_MAX_DELAY = 7

# In days, how many day in the past to look at when initializing database
# for regular paper consumer
CONSUMER_INIT_PAST = 365
# for pubstream
CONSUMER_PUBPEER_INIT_PAST = 365

CONSUMER_PUBMED_EMAIL = env('CONSUMER_PUBMED_EMAIL')
CONSUMER_ELSEVIER_API_KEY = env('CONSUMER_ELSEVIER_API_KEY')
CONSUMER_ELSEVIER_API_KEY2 = env('CONSUMER_ELSEVIER_API_KEY2')
CONSUMER_SPRINGER_KEY = env('SPRINGER_KEY')
CONSUMER_PUBPEER_API_KEY = env('PUBPEER_KEY')
CONSUMER_PUBPEER_USER_EMAIL = 'contact@pubpeer.com'


# NLP APP
# ------------------------------------------------------------------------------
NLP_CHUNK_SIZE = 10000
NLP_DATA_PATH = str((ROOT_DIR-1).path('data', 'nlp', 'data'))
NLP_MODELS_PATH = str((ROOT_DIR-1).path('data', 'nlp', 'mods'))
NLP_PE_PATH = str((ROOT_DIR-1).path('data', 'nlp', 'pe'))
NLP_TE_PATH = str((ROOT_DIR-1).path('data', 'nlp', 'te'))
NLP_NLTK_DATA_PATH = str(APPS_DIR.path('nlp', 'nltk_data'))
NLP_MAX_VECTOR_SIZE = 300
NLP_MAX_KNN_NEIGHBORS = 10

# Time in days for recomputing neighbors is accessed
NLP_NEIGHBORS_REFRESH_TIME_LAPSE = 7


# LIBRARY APP
# ------------------------------------------------------------------------------

LIBRARY_ITEMS_PER_PAGE = 100
LIBRARY_NUMBER_OF_NEIGHBORS = 5
LIBRARY_NUMBER_OF_THREADS_NEIGHBORS = 5
LIBRARY_DEFAULT_NEIGHBORS_TIMESPAN = 60


# THREAD APP
# ------------------------------------------------------------------------------

MAX_MEMBERS_NOT_JOINED = 3
THREADS_NUMBER_OF_NEIGHBORS = 5
THREADS_DEFAULT_NEIGHBORS_TIMESPAN = 60


# FEED APP
# ------------------------------------------------------------------------------
# Feeds default controls states
FEEDS_DEFAULT_TIME_SPAN = 30

HIDE_CLUSTER_ICON = True

FEED_TIME_SPAN_DEFAULT = 30

FEED_N_FIRST_PAPERS_ONLY = 1000

FEED_STREAM_SCORE_THRESHOLD_DEFAULT = 0.3
FEED_TREND_SCORE_THRESHOLD_DEFAULT = 0.1
FEED_THREADFEED_SCORE_THRESHOLD_DEFAULT = 0.3

# ALTMETRIC APP
# ------------------------------------------------------------------------------
ALTMETRIC_API_KEY = env('ALTMETRIC_API_KEY')
# slightly less than each second in a a day
ALTMETRIC_MAX_PAPERS_PER_PERIOD = 20 * 3600


# AVATAR
# ------------------------------------------------------------------------------
AVATAR_AUTO_GENERATE_SIZES = (80,)
AVATAR_CACHE_ENABLED = True
AVATAR_DEFAULT_URL = 'https://s3-us-west-2.amazonaws.com/etalia-production-static/avatar.jpg'
AVATAR_DEFAULT_SIZE = 80
AVATAR_EXPOSE_USERNAMES = False
AVATAR_GRAVATAR_DEFAULT = None
AVATAR_GRAVATAR_FORCEDEFAULT = False
AVATAR_GRAVATAR_FIELD = 'email'
AVATAR_MAX_SIZE = 1024 * 1024
AVATAR_STORAGE_DIR = 'photos'
AVATAR_ADD_TEMPLATE = 'avatar/add.html'
AVATAR_CHANGE_TEMPLATE = 'avatar/change.html'
AVATAR_DELETE_TEMPLATE = 'avatar/confirm_delete.html'


# POPOVERS
# ------------------------------------------------------------------------------
# Display All highest priority if POPOVERS_DISPLAY_HIGHEST_PRIORITY is True
# plus POPOVERS_DISPLAY_NEW_ANCHORED NEW modal type
# plus POPOVERS_DISPLAY_NEW_ANCHORED NEW anchored type

POPOVERS_DISPLAY_HIGHEST_PRIORITY = True
POPOVERS_DISPLAY_NEW_ANCHORED = 1
POPOVERS_DISPLAY_NEW_MODAL = 1
POPOVERS_DISPLAY_REFRESH_PERIOD = 3600  # in seconds


# EASY TIMEZONE
# ------------------------------------------------------------------------------
GEOIP_DATABASE = os.path.join(str(ROOT_DIR), 'data/GeoIP/GeoLiteCity.dat')
GEOIPV6_DATABASE = os.path.join(str(ROOT_DIR), 'data/GeoIP/GeoLiteCityv6.dat')

# REST FRAMEWORK
# ------------------------------------------------------------------------------
REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': (
        'rest_framework.permissions.IsAuthenticated',
    ),
    'DEFAULT_METADATA_CLASS': 'rest_framework.metadata.SimpleMetadata',
    'DEFAULT_PAGINATION_CLASS': 'rest_framework.pagination.LimitOffsetPagination',
    'PAGE_SIZE': 10,
    'URL_FIELD_NAME': 'link',
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'rest_framework.authentication.SessionAuthentication',
        'rest_framework.authentication.TokenAuthentication',
    ),
    'DEFAULT_THROTTLE_CLASSES': (
        'rest_framework.throttling.AnonRateThrottle',
        'rest_framework.throttling.UserRateThrottle'
    ),
    'DEFAULT_THROTTLE_RATES': {
        'anon': '1000/day',
        'user': '100000/day'
    }
}


# CACHE
# ------------------------------------------------------------------------------
CACHE_MIDDLEWARE_SECONDS = 60 * 5
CACHE_MIDDLEWARE_ALIAS = 'default'
CACHE_MIDDLEWARE_KEY_PREFIX = ''


# DJANGO TOOLBAR
# ------------------------------------------------------------------------------
DEBUG_TOOLBAR = False


# LOGGING CONFIGURATION
# ------------------------------------------------------------------------------
LOG_DIR = (ROOT_DIR-1).path('log')
# creating log directory if does not exist
if not os.path.exists(str(LOG_DIR)):
    os.mkdir(str(LOG_DIR))

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
            'filename': str(LOG_DIR.path('etalia.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'populate': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('populate.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'nlp': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('nlp.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'feeds': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('feeds.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'users': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('users.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'consumers': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('consumers.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'celery': {
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': str(LOG_DIR.path('celery.log')),
            'formatter': 'verbose',
            'maxBytes': 1024 * 1024 * 5,  # 5 mb
        },
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler',
            'include_html': True,
        },
    },
    'loggers': {
        'django': {
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': False,
        },
        'etalia': {
            'handlers': ['console', 'file'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
        },
        'etalia.populate': {
            'handlers': ['console', 'populate'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'etalia.consumers': {
            'handlers': ['console', 'consumers', 'mail_admins'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'etalia.nlp': {
            'handlers': ['console', 'nlp', 'mail_admins'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'etalia.feeds': {
            'handlers': ['console', 'feeds', 'mail_admins'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
        'etalia.users': {
            'handlers': ['console', 'users', 'mail_admins'],
            'level': env.str('DJANGO_LOG_LEVEL', 'INFO'),
            'propagate': False,
        },
    }
}
