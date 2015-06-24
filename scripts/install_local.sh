#!/usr/bin/env bash

# install local requirement
pip install -r ../requirements/local.txt

## Set environment variable for current virtual environment
#./routines/init_env_vars.sh
#
## Create Postgres database
#./routines/create_db.sh
#
## Init migrations and migrate
#./routines/init_migrations.sh
#
## Populate database
#./routines/init_populate.sh

# you should now be able to fire up the server locally
# > ./manage.py runserver --settings=config.settings.local
#
# However, for asynchroneous tasks to run you need to fire up celery and
# rabbitmq-server up front
# > rabbitmq-server &
# > celery -A config worker --loglevel=info