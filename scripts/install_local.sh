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
: "${DJANGO_SETTINGS_MODULE?Need to set DJANGO_SETTINGS_MODULE}"
: "${DJANGO_LOG_LEVEL?Need to set DJANGO_LOG_LEVEL}"
: "${DJANGO_DEBUG?Need to set DJANGO_DEBUG}"
: "${CONSUMER_ELSEVIER_API_KEY?Need to set CONSUMER_ELSEVIER_API_KEY}"
: "${CONSUMER_PUBMED_EMAIL?Need to set CONSUMER_PUBMED_EMAIL}"
: "${SOCIAL_AUTH_MENDELEY_KEY?Need to set SOCIAL_AUTH_MENDELEY_KEY}"
: "${SOCIAL_AUTH_MENDELEY_SECRET?Need to set SOCIAL_AUTH_MENDELEY_SECRET}"
: "${SOCIAL_AUTH_ZOTERO_KEY?Need to set SOCIAL_AUTH_ZOTERO_KEY}"
: "${SOCIAL_AUTH_ZOTERO_SECRET?Need to set SOCIAL_AUTH_ZOTERO_SECRET}"
: "${DJANGO_AWS_REGION?Need to set DJANGO_AWS_REGION}"
: "${DJANGO_AWS_ACCESS_KEY_ID?Need to set DJANGO_AWS_ACCESS_KEY_ID}"
: "${DJANGO_AWS_SECRET_ACCESS_KEY?Need to set DJANGO_AWS_SECRET_ACCESS_KEY}"
: "${DJANGO_AWS_STORAGE_BUCKET_NAME?Need to set DJANGO_AWS_STORAGE_BUCKET_NAME}"
: "${DJANGO_SETTINGS_MODULE?Need to set DJANGO_SETTINGS_MODULE}"
: "${DISQUS_WEBSITE_SHORTNAME?Need to set DISQUS_WEBSITE_SHORTNAME}"
: "${DISQUS_PUBLIC_KEY?Need to set DISQUS_PUBLIC_KEY}"
: "${DISQUS_SECRET_KEY?Need to set DISQUS_SECRET_KEY}"
: "${ALTMETRIC_API_KEY?Need to set ALTMETRIC_API_KEY}"


# MAKE FOLDER
mkdir ../logs

# MAKE MIGRATIONS AND MIGRATE
../manage.py makemigrations
../manage.py migrate

# POPULATE DATABASE
# populate library with some test data
../manage.py populate publisher all
../manage.py populate journal thomson_local
../manage.py populate journal pubmed_local
../manage.py populate journal arxiv_local
# populate consumers
../manage.py populate consumer pubmed --name pubmed_all --local
../manage.py populate consumer arxiv --name arxiv_all
../manage.py populate consumer elsevier --name elsevier_all

# fetch some new papers
./manage.py shell < init.py

# RUN TEST
py.test

