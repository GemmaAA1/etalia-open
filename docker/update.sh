#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OPTIONS="-i --rm \
    --volumes-from etalia_web \
    --network etalia-network \
    --link etalia_db:etalia_db \
    --link etalia_rabbit:etalia_rabbit \
    --link etalia_redis:etalia_redis \
    --env-file $DIR/.envs \
    etalia/python-dev"

docker-compose -f ./docker/docker-compose.yml up -d
docker run $OPTIONS python setup/manager.py --update