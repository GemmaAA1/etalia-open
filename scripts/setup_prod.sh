#!/usr/bin/env bash

#!/usr/bin/env bash
export DJANGO_SETTINGS_MODULE='config.settings.production'
export DJANGO_LOG_LEVEL='DEBUG'
export DJANGO_DEBUG=True
export DJANGO_SECRET_KEY='_kjjr3)+tdlfbcw7uu&oue+*50+hbv9gsd-yx35^*%n$5ugp-s'
export CONSUMER_ELSEVIER_API_KEY='293e288c325d7765b7c22f5195175351'
export CONSUMER_PUBMED_EMAIL='nicolas.pannetier@gmail.com'
export SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_KEY='1678'
export SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_SECRET='caOrLU0DqOUC4wdD'
export SOCIAL_AUTH_CUSTOM_ZOTERO_KEY='a7ecbff3d0bbe59abc4b'
export SOCIAL_AUTH_CUSTOM_ZOTERO_SECRET='c5d0c178d9196e62bdbf'
export DJANGO_AWS_ACCESS_KEY_ID='AKIAJP4QVWCJZTBCDW7A'
export DJANGO_AWS_SECRET_ACCESS_KEY='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'
export DJANGO_AWS_STORAGE_BUCKET_NAME='paperstreamstatic'
export DISQUS_WEBSITE_SHORTNAME='paperstream'
export DISQUS_PUBLIC_KEY='w2W0iBEJwGE49PjupwQxDnfzC9ayliEvctiGwbmVb63uHIXNTZLgreJDNRvvBOap'
export DISQUS_SECRET_KEY='eMWsm6qeNkDHzdvLViScWPldyDVnmvAz4U79YjsCelOu58XnRPelrUimqTrhGrRw'
export AWS_RDS_DB_NAME='paperstream_staging'
export AWS_RDS_USERNAME='npannetier'
export AWS_RDS_PASSWORD='Octopusisdb'
export AWS_RDS_HOSTNAME='paperstream-db-staging.c61ni52g1qt1.us-west-2.rds.amazonaws.com'
export AWS_RDS_PORT='5432'

# CREATE Postgres DATABASE
psql postgres -c "CREATE DATABASE paperstream"

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

