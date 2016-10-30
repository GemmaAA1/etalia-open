# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### Setting things up for development:

Etalia stack: Django, Postgres, RabbitMQ, Redis and Celery. 
To setup etalia locally you need [Docker](https://www.docker.com/).   

1. Clone Etalia repository:
        
        git clone https://[your_bitbucket_username]@bitbucket.org/NPann/etalia.git

2. Run install.sh

        $ ./install.sh

3. Launch containers:        

        $ docker-compose up

4. visit 127.0.0.1:8000


### Updating Etalia database and engines:

To update etalia database and engines run:

        $ docker-compose run full ./docker/update.sh

### Frontend dev ###

During development, use ```gulp``` to copy librairies from *nodes_modules/* and *bower_components/* to *static/js/lib/*.

For production, use ```gulp prod``` to build compiled assets. 

New librairies must be configured in *static/js/*

### Etalia team ###

* Etienne Dauvergne
* Nicolas Pannetier
* Norbert Schuff
