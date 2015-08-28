from .common import *

DEBUG = env.bool('DJANGO_DEBUG', default=False)

# Parse database configuration from $DATABASE_URL
import dj_database_url
DATABASES['default'] =  dj_database_url.config()

# Honor the 'X-Forwarded-Proto' header for request.is_secure()
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

# Allow all host headers
ALLOWED_HOSTS = ['*']

# Static asset configuration
import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_ROOT = 'staticfiles'
STATIC_URL = '/static/'

STATICFILES_DIRS = (
    # os.path.join(BASE_DIR, 'static'),
)


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

