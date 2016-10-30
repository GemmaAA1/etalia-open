#!/bin/sh

docker network create etalia-network
docker volume create --name=etalia-data-volume

docker-compose -f ./docker/docker-compose.yml up -d

./docker/migrate.sh
./docker/load.sh

docker restart etalia_worker

./docker/update.sh
