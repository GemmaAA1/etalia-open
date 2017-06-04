#!/bin/sh

docker network create etalia-network
docker volume create --name=etalia-data-volume

docker-compose -f ./docker/docker-compose.yml up -d

./migrate.sh
./load.sh

docker restart etalia_worker

./update.sh