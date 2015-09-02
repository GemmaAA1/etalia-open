# -*- coding: utf-8 -*-
from .common import *

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['*']

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': env.str('RDS_DB_NAME'),
        'USER': env.str('RDS_USERNAME'),
        'PASSWORD': env.str('RDS_PASSWORD'),
        'HOST': env.str('RDS_HOSTNAME'),
        'PORT': env.str('RDS_PORT'),
    }
}