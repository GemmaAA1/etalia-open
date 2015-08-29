# -*- coding: utf-8 -*-

from .heroku_local import *

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['shrouded-peak-3632.herokuapp.com',
                 'shrouded-peak.paperstream.io']

# Parse database configuration from $DATABASE_URL
import dj_database_url
DATABASES['default'] = dj_database_url.config()
