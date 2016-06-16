#!/usr/bin/env bash

# This is script must be run from the virtual env you have defined.
# -> python3.4, pip3, psql (Postgres), rabbitmq are considered installed
#
# After running that script you should be able to fire up the server locally
# > ./manage.py runserver --settings=config.settings.local
#
# For asynchroneous tasks to run you need to fire up celery and
# rabbitmq-server up front
# > rabbitmq-server &
# > celery -A config worker --loglevel=info
#
# Also modified 'USER' in DATABASES definition in config/common.py to match
# you postgres username.

# REQUIREMENTS
pip3 install -r ../requirements/local.txt

# ENVIRONMENT VARIABLE IN VIRTUAL ENV
# Check if virtual env is active
if [ -z "$VIRTUAL_ENV" ]; then
    echo 'VIRTUAL_ENV not set'
    exit 1
fi

# MAKE MIGRATIONS AND MIGRATE
../manage.py makemigrations
../manage.py migrate

# POPULATE DATABASE
../manage.py shell < routines/update.py --init-production

# RUN TEST
#py.test

