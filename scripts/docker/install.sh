#!/bin/sh

docker network create etalia-network
docker volume create --name=etalia-data-volume

DOCKER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../docker" && pwd )"
docker-compose -f ${DOCKER_DIR}/docker-compose.yml up -d

./migrate.sh
./load.sh

docker restart etalia_worker

./update.sh