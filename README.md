# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### Setting things up for development:

Etalia stack: Django, Postgres, RabbitMQ, Redis and Celery. 
To setup etalia locally you need [Docker](https://www.docker.com/).   

1. Clone Etalia repository:
        
        git clone https://[your_bitbucket_username]@bitbucket.org/NPann/etalia.git


2. Run setup.sh (this will build, run and initialize the Etalia stack):

        $ ./setup.sh

3. Visit _127.0.0.1:8000_

4. Read more in [docker doc](./docs/docker.md).

NB: If you get S3Forbidden errors, likely your docker images clock went off sync. Run ./docker/sync_clock.sh to sync
 them back. If you are still having a S3Forbidden error, it is something else...


### Frontend dev ###

During development, use ```gulp``` to copy libraries from *nodes_modules/* and *bower_components/* to *static/js/lib/*.

For production, use ```gulp prod``` to build compiled assets. 

New librairies must be configured in *static/js/*

### Etalia team ###

* Etienne Dauvergne
* Nicolas Pannetier
* Norbert Schuff
