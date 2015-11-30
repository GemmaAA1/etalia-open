# README #

PubStream is an webapp that provides scientific publication readers with relevant
 recent scientific publications.

### How do I get set up locally? ###

PubStream stack is based on:
Django
PostgreSQL
RabbitMQ
Celery
 
Pre-requisites are: 
* pip (the python package management system)
* virtualenv
* PostgreSQL
* RabbitMQ
* Some default environment variables (ask for it)

Update configuration file with local database
* Create a local Postgres database (e.g. "CREATE pubstream_dev;" from psql on mac)
* Update paperstream/config/settings/common.py with your username and database name
   ex: 
   DATABASES = {
      'default': {
          ...
          'NAME': 'localdatabasename',
          'USER': 'yourpostgresusername',
          'PASSWORD': '',
          ...
      }
   }
* Run install script  

```
#!bash
$ ./scripts/setup_local.sh
```

* Fire-up the server:

```
#!bash
$ ./manage.py runserver
```

* Visit 127.0.0.1:8000

### Who do I talk to? ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)


