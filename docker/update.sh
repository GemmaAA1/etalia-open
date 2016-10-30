#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OPTIONS="-i --rm \
    --volumes-from etalia_web \
    --network etalia-network \
    --link etalia_db:db \
    --link etalia_rabbit:rabbit \
    --link etalia_redis:redis \
    --env-file $DIR/.envs \
    etalia/python-dev"

docker run $OPTIONS python setup/manager.py --update
