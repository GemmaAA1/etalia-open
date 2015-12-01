# README #

PubStream is an webapp that provides scientific publication readers with relevant
 recent scientific publications.

### How do I get set up locally? ###

PubStream stack is based on: Django, PostgreSQL, RabbitMQ, Celery

Pre-requisites are: pip (the python package management system), virtualenv, PostgreSQL, RabbitMQ.

Once you get the pre-requisites installed, follow these steps to get started:

1. Create you virtual environment and export environment variable from the development configuration file (Request one if you don't have one already)
2. Clone PubStream: 
```
#!bash

git clone https://NPann@bitbucket.org/NPann/paperstream.git
```

3. Create a local Postgres database (e.g. "CREATE DATABASE pubstream_dev;" from psql on mac)
4. Update database configuration in paperstream/config/settings/common.py with your username and database name
   ex: 
   DATABASES = {
      'default': {
          ...
          'NAME': 'database_name',
          'USER': 'your_username',
          'PASSWORD': '',
          ...
      }
   }
   
3. Run install_local.sh script from directory paperstream/scripts

```
#!bash

$ ./install_local.sh
```

* Fire-up the server:

```
#!bash

$ ./manage.py runserver
```

* Visit 127.0.0.1:8000

### Who do I talk to? ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)


