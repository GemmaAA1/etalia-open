#!/bin/sh

DOCKER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../docker" && pwd )"

OPTIONS="-i --rm \
    --volumes-from etalia_web \
    --network etalia-network \
    --link etalia_db:db \
    --env-file $DOCKER_DIR/.envs \
    etalia/dev"

docker run $OPTIONS python manage.py makemigrations
docker run $OPTIONS python manage.py migrate
