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
# populate library
../manage.py populate publisher all
../manage.py populate journal thomson
../manage.py populate journal pubmed
../manage.py populate journal arxiv
# populate consumers
../manage.py populate consumer pubmed --name pubmed_all
../manage.py populate consumer arxiv --name arxiv_all
../manage.py populate consumer elsevier --name elsevier_all

# RUN TEST
py.test

