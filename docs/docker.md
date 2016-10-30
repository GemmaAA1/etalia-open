Docker
======

First of all, create a volume to store the database:

    $ docker volume create --name=etalia_db_vol

To build and run the containers, go to _./docker_ directory and run docker compose up:

    $ cd [project-path]/docker
    $ docker-compose up -d

To remove the containers, go to _./docker_ directory and run docker compose down:

    $ cd [project-path]/docker
    $ docker-compose down


https://docs.docker.com/engine/reference/commandline/
https://docs.docker.com/compose/compose-file/

https://github.com/andrecp/django-tutorial-docker-nginx-postgres
https://realpython.com/blog/python/django-development-with-docker-compose-and-machine/
http://www.fullstackpython.com/docker.html
https://www.digitalocean.com/community/tutorials/docker-explained-how-to-containerize-python-web-applications
