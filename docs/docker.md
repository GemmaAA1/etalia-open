Docker
======

First of all, create a volume to store the database, and a network:

    $ docker volume create --name=etalia-data-volume
    $ docker network create etalia-network

To build the containers, go to _./docker_ directory and run docker compose build:

    $ cd [project-path]/docker
    $ docker-compose build
                            
To launch the containers run:

    $ cd [project-path]/docker
    $ docker-compose up -d

To run migrations of etalia database:

    $ ./docker/migrate.sh

To load the data fixtures:

    $ ./docker/load.sh

To update etalia data and engines

    $ ./docker/update.sh
                    
To stop the containers, go to _./docker_ directory and run docker compose stop:

    $ cd [project-path]/docker
    $ docker-compose stop

To remove the containers, go to _./docker_ directory and run docker compose down:

    $ cd [project-path]/docker
    $ docker-compose down -v --remove-orphans


https://docs.docker.com/engine/reference/commandline/
https://docs.docker.com/compose/compose-file/

https://github.com/andrecp/django-tutorial-docker-nginx-postgres
https://realpython.com/blog/python/django-development-with-docker-compose-and-machine/
http://www.fullstackpython.com/docker.html
https://www.digitalocean.com/community/tutorials/docker-explained-how-to-containerize-python-web-applications
