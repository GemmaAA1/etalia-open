#!/bin/sh

DOCKER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../docker" && pwd )"

OPTIONS="-i --rm \
    -v $DOCKER_DIR/../.:/code
    --network etalia-network \
    --link etalia_db:etalia_db \
    --link etalia_rabbit:etalia_rabbit \
    --link etalia_redis:etalia_redis \
    --env-file $DOCKER_DIR/.envs \
    etalia/dev"

docker run $OPTIONS python setup/manager.py --update