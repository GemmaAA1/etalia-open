from .dev import *

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'test',
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
