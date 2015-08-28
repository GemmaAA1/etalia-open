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

# install local requirement
pip3 install -r ../requirements/local.txt

# Set environment variable for current virtual environment
./routines/init_env_vars.sh

# Create local files
./routines/init_files.sh

# Create Postgres database
./routines/create_db.sh

# Init migrations and migrate
./routines/init_migrations.sh

# Populate database
./routines/init_populate.sh

# download nltk_data
python ./routines/init_nltk.py

# test
py.test

