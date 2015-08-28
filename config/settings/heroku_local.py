from .heroku import *

# Restore Dev DATABASES setup
from .dev import DATABASES as DB_DEV
DATABASES['default'] = DB_DEV['default']

ALLOWED_HOSTS = ALLOWED_HOSTS + ['localhost', '0.0.0.0']
