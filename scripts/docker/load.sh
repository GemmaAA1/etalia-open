#!/bin/sh

DOCKER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../docker" && pwd )"

OPTIONS="-i --rm \
    -v $DOCKER_DIR/../.:/code
    --network etalia-network \
    --link etalia_db:db \
    --env-file $DOCKER_DIR/.envs \
    etalia/dev"



docker run $OPTIONS python control/manager.py --load
