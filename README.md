# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### Setting things up for development:

For development, Etalia requires:

* Python 3.4+
* pip (the python package manager)
* virtualenv
* PostgreSQL
* RabbitMQ

Once you get the requirements installed, follow these steps to get started:

1. Create you virtual environment with Python 3.4 onboard. Export the environment variables for Etalia (Request them if you don't have them). 

2. Clone Etalia repository:
        
        git clone https://[your_bitbucket_username]@bitbucket.org/NPann/etalia.git

3. Create a local Postgres database (e.g. `CREATE DATABASE etalia;` from psql on mac)

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
