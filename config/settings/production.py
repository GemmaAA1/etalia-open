from .base import *

DEBUG = env.bool('DJANGO_DEBUG', default=False)
SECRET_KEY = env('PAP_SECRET_KEY')

# Parse database configuration from $DATABASE_URL
import dj_database_url
DATABASES['default'] = dj_database_url.config()

# Honor the 'X-Forwarded-Proto' header for request.is_secure()
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')

# Allow all host headers
ALLOWED_HOSTS = ['*']

# Static asset configuration
import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_ROOT = 'staticfiles'
STATIC_URL = '/static/'
STATICFILES_DIRS = []

# Amazon
# S3
INSTALLED_APPS += ('storages',)
AWS_STORAGE_BUCKET_NAME = "paperstreamstatic"
STATICFILES_STORAGE = 'storages.backends.s3boto.S3BotoStorage'
S3_URL = 'http://%s.s3.amazonaws.com/' % AWS_STORAGE_BUCKET_NAME
STATIC_URL = S3_URL

# RDS
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'paperstream_default',
        'USER': 'nicolaspannetier',
        'PASSWORD': 'db_Oct0pu$',
        'HOST': 'paperstreamdb.cxnh00eqsn0h.us-west-2.rds.amazonaws.com:5432',
        'PORT': '5432',
    }
}