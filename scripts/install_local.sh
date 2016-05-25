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
pip3 install -r ../requirements/development.txt
# upgrade numpy (error between 1.8.2 and gensim (does not seem the case on ubuntu...))
pip install --upgrade numpy

# VIRTUAL ENV
# Check if virtual env is active
if [ -z "$VIRTUAL_ENV" ]; then
    echo 'VIRTUAL_ENV not set'
    exit 1
fi

# ENVIRONMENT VARIABLES
# Check if required env variable are defined
: "${DJANGO_SETTINGS_MODULE?Error: Need to set env var DJANGO_SETTINGS_MODULE}"
: "${DJANGO_LOG_LEVEL?Error: Need to set env var DJANGO_LOG_LEVEL}"
: "${DJANGO_DEBUG?Error: Need to set env var DJANGO_DEBUG}"
: "${CONSUMER_ELSEVIER_API_KEY?Error: Need to set env var CONSUMER_ELSEVIER_API_KEY}"
: "${CONSUMER_PUBMED_EMAIL?Error: Need to set env var CONSUMER_PUBMED_EMAIL}"
: "${SOCIAL_AUTH_MENDELEY_KEY?Error: Need to set env var SOCIAL_AUTH_MENDELEY_KEY}"
: "${SOCIAL_AUTH_MENDELEY_SECRET?Error: Need to set env var SOCIAL_AUTH_MENDELEY_SECRET}"
: "${SOCIAL_AUTH_ZOTERO_KEY?Error: Need to set env var SOCIAL_AUTH_ZOTERO_KEY}"
: "${SOCIAL_AUTH_ZOTERO_SECRET?Error: Need to set env var SOCIAL_AUTH_ZOTERO_SECRET}"
: "${AWS_REGION?Error: Need to set env var AWS_REGION}"
: "${AWS_ACCESS_KEY_ID?Error: Need to set env var AWS_ACCESS_KEY_ID}"
: "${AWS_SECRET_ACCESS_KEY?Error: Need to set env var AWS_SECRET_ACCESS_KEY}"
: "${AWS_STORAGE_BUCKET_NAME?Error: Need to set env var AWS_STORAGE_BUCKET_NAME}"
: "${DJANGO_SETTINGS_MODULE?Error: Need to set env var DJANGO_SETTINGS_MODULE}"
: "${DISQUS_WEBSITE_SHORTNAME?Error: Need to set env var DISQUS_WEBSITE_SHORTNAME}"
: "${DISQUS_PUBLIC_KEY?Error: Need to set env var DISQUS_PUBLIC_KEY}"
: "${DISQUS_SECRET_KEY?Error: Need to set env var DISQUS_SECRET_KEY}"
: "${ALTMETRIC_API_KEY?Error: Need to set env var ALTMETRIC_API_KEY}"

# init database with papers, nlp models and altmetric data
../manage.py shell < routines/update.py --init

# RUN TEST
# py.test

