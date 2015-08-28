import os
from __future__ import absolute_import, unicode_literals
from boto.s3.connection import OrdinaryCallingFormat
from .common import *

DEBUG = env.bool('DJANGO_DEBUG', default=False)

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['localhost', 'paperstream.herokuapp.com']

# Parse database configuration from $DATABASE_URL
import dj_database_url
DATABASES['default'] = dj_database_url.config()

# Honor the 'X-Forwarded-Proto' header for request.is_secure()
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

# Static asset configuration
STATIC_ROOT = 'staticfiles'
STATIC_URL = '/static/'

STATICFILES_DIRS = (
    str(APPS_DIR.path('static')),
)

STATICFILES_STORAGE = 'whitenoise.django.GzipManifestStaticFilesStorage'

# # Static asset cocnfiguration
# STATIC_ROOT = 'staticfiles'
# STATICFILES_DIRS = []
#
# # Amazon
# # S3
# INSTALLED_APPS += ('storages',)
# AWS_STORAGE_BUCKET_NAME = "paperstreamstatic"
# STATICFILES_STORAGE = 'storages.backends.s3boto.S3BotoStorage'
# S3_URL = 'http://%s.s3.amazonaws.com/' % AWS_STORAGE_BUCKET_NAME
# STATIC_URL = S3_URL


# # STORAGE CONFIGURATION
# # ------------------------------------------------------------------------------
# # Uploaded Media Files
# # ------------------------
# # See: http://django-storages.readthedocs.org/en/latest/index.html
# INSTALLED_APPS += (
#     'storages',
# )
# DEFAULT_FILE_STORAGE = 'storages.backends.s3boto.S3BotoStorage'
#
# AWS_ACCESS_KEY_ID = env('DJANGO_AWS_ACCESS_KEY_ID')
# AWS_SECRET_ACCESS_KEY = env('DJANGO_AWS_SECRET_ACCESS_KEY')
# AWS_STORAGE_BUCKET_NAME = env('DJANGO_AWS_STORAGE_BUCKET_NAME')
# AWS_AUTO_CREATE_BUCKET = True
# AWS_QUERYSTRING_AUTH = False
# AWS_S3_CALLING_FORMAT = OrdinaryCallingFormat()
#
# # AWS cache settings, don't change unless you know what you're doing:
# AWS_EXPIRY = 60 * 60 * 24 * 7
#
# # TODO See: https://github.com/jschneier/django-storages/issues/47
# # Revert the following and use str after the above-mentioned bug is fixed in
# # either django-storage-redux or boto
# AWS_HEADERS = {
#     'Cache-Control': six.b('max-age=%d, s-maxage=%d, must-revalidate' % (
#         AWS_EXPIRY, AWS_EXPIRY))
# }
#
# # URL that handles the media served from MEDIA_ROOT, used for managing
# # stored files.
# MEDIA_URL = 'https://s3.amazonaws.com/%s/' % AWS_STORAGE_BUCKET_NAME
#
# # Static Assests
# # ------------------------
# {% if cookiecutter.use_whitenoise == 'y' -%}
# STATICFILES_STORAGE = 'whitenoise.django.GzipManifestStaticFilesStorage'
# {% else %}
# STATICFILES_STORAGE = DEFAULT_FILE_STORAGE
# STATIC_URL = MEDIA_URL
#
# # See: https://github.com/antonagestam/collectfast
# # For Django 1.7+, 'collectfast' should come before
# # 'django.contrib.staticfiles'
# AWS_PRELOAD_METADATA = True
# INSTALLED_APPS = ('collectfast', ) + INSTALLED_APPS
# {%- endif %}
