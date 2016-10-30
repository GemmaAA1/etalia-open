#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OPTIONS="-i --rm \
    --volumes-from etalia_web \
    --network etalia-network \
    --link etalia_db:db \
    --env-file $DIR/.envs \
    etalia/python-dev"

docker run $OPTIONS python manage.py makemigrations
docker run $OPTIONS python manage.py migrate
