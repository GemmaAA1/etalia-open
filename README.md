# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### Setting things up for development:

Etalia stack: Django, Postgres, RabbitMQ, Redis and Celery. 
To setup etalia locally you need [Docker](https://www.docker.com/).   

1. Clone Etalia repository:
        
        git clone https://[your_bitbucket_username]@bitbucket.org/NPann/etalia.git

2. Build the docker image (must be done each time the python dependencies changes):

        $ ./docker/build.sh

3. Run install.sh (this will build, run and initialize the Etalia stack):

        $ ./install.sh

4. Visit _127.0.0.1:8000_

5. Read more in [docker doc](./docs/docker.md).

### Frontend dev ###

During development, use ```gulp``` to copy librairies from *nodes_modules/* and *bower_components/* to *static/js/lib/*.

For production, use ```gulp prod``` to build compiled assets. 

New librairies must be configured in *static/js/*

### Etalia team ###

* Etienne Dauvergne
* Nicolas Pannetier
* Norbert Schuff
