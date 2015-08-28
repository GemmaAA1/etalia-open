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
#export environment variable to virtulalenv postactivate
file=$VIRTUAL_ENV/bin/postactivate
rm $file
cat >$file <<EOL
#!/usr/bin/env bash
export DJANGO_SETTINGS_MODULE='config.settings.dev'
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
EOL
source $file

## cleaning up symmetrically
file=$VIRTUAL_ENV/bin/predeactivate
rm $file
cat >$file <<EOL
#!/usr/bin/env bash
unset DJANGO_SETTINGS_MODULE='config.settings.dev'
unset DJANGO_LOG_LEVEL='DEBUG'
unset DJANGO_DEBUG=True
unset DJANGO_SECRET_KEY='_kjjr3)+tdlfbcw7uu&oue+*50+hbv9gsd-yx35^*%n$5ugp-s'
unset CONSUMER_ELSEVIER_API_KEY='293e288c325d7765b7c22f5195175351'
unset CONSUMER_PUBMED_EMAIL='nicolas.pannetier@gmail.com'
unset SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_KEY='1678'
unset SOCIAL_AUTH_CUSTOM_MENDELEY_OAUTH2_SECRET='caOrLU0DqOUC4wdD'
unset SOCIAL_AUTH_CUSTOM_ZOTERO_KEY='a7ecbff3d0bbe59abc4b'
unset SOCIAL_AUTH_CUSTOM_ZOTERO_SECRET='c5d0c178d9196e62bdbf'
unset DJANGO_AWS_ACCESS_KEY_ID='AKIAJP4QVWCJZTBCDW7A'
unset DJANGO_AWS_SECRET_ACCESS_KEY='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'
unset DJANGO_AWS_STORAGE_BUCKET_NAME='paperstreamstatic'
unset DISQUS_WEBSITE_SHORTNAME='paperstream'
unset DISQUS_PUBLIC_KEY='w2W0iBEJwGE49PjupwQxDnfzC9ayliEvctiGwbmVb63uHIXNTZLgreJDNRvvBOap'
unset DISQUS_SECRET_KEY='eMWsm6qeNkDHzdvLViScWPldyDVnmvAz4U79YjsCelOu58XnRPelrUimqTrhGrRw'
EOL

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

