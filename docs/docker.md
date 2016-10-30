Docker
======

First of all, create a volume to store the database:

    $ docker volume create --name=etalia-data-volume

To build the containers, go to _./docker_ directory and run docker compose build:

    $ cd [project-path]/docker
    $ docker-compose build

To run migrations of etalia database:

    $ docker-compose run simple ./docker/migrate.sh

To update etalia data and engines

    $ docker-compose run full ./docker/update.sh

To launch containers run:

    $ cd [project-path]/docker
    $ docker-compose up

To stop the containers, go to _./docker_ directory and run docker compose stop:

    $ cd [project-path]/docker
    $ docker-compose stop

To remove the containers, go to _./docker_ directory and run docker compose down:

    $ cd [project-path]/docker
    $ docker-compose down


https://docs.docker.com/engine/reference/commandline/
https://docs.docker.com/compose/compose-file/

https://github.com/andrecp/django-tutorial-docker-nginx-postgres
https://realpython.com/blog/python/django-development-with-docker-compose-and-machine/
http://www.fullstackpython.com/docker.html
https://www.digitalocean.com/community/tutorials/docker-explained-how-to-containerize-python-web-applications
