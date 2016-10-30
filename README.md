# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### Setting things up for development:

Etalia stack is composed of: Django, Postgres, RabbitMQ, Redis and Celery workers. 
To setup etalia locally you need [Docker](https://www.docker.com/).   

1. Clone Etalia repository:
        
        git clone https://[your_bitbucket_username]@bitbucket.org/NPann/etalia.git

2. Create a etalia volume:

        $ docker volume create --name=etalia_db_vol

3. Go to _./docker_ directory and run docker-compose up:

        $ cd [project-path]/docker
        $ docker-compose up -d

4. Copy etalia/config/settings/common.py.dist to etalia/config/settings/common.py

5. Update database settings in etalia/config/settings/common.py with your username and database name. Example:
    
        DATABASES = {
            'default': {
                ...
                'NAME': 'database_name',
                'USER': 'your_username',
                'PASSWORD': '',
                ...
            }
        }

6. Go to scripts directory and run the install_local script:

        cd scripts/
        ./install_local.sh

### Launching the app:

`
./manage.py runserver
`

and visit 127.0.0.1:8000

### Updating Etalia:

An update script is provided to update database and related objects.
See help menu for usage:
`
./scripts/routines/update.py -h
`

### Frontend dev ###

During development, use ```gulp``` to copy librairies from *nodes_modules/* and *bower_components/* to *static/js/lib/*.

For production, use ```gulp prod``` to build compiled assets. 

New librairies must be configured in *static/js/*

### Etalia team ###

* Etienne Dauvergne
* Nicolas Pannetier (nicolas.pannetier@gmail.com)
* Norbert Schuff
