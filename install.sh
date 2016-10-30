#!/bin/sh
docker volume create --name=etalia-data-volume
cd docker
docker-compose build
docker-compose run simple ./docker/migrate.sh
docker-compose run simple ./docker/load.sh
docker restart etalia_worker
docker-compose run full ./docker/update.sh