# README #

PubStream is an webapp that provides scientific publication readers with relevant
 recent scientific publications.

### How do I get set up? ###

PubStream stack is based on:
Django
PostgreSQL
RabbitMQ
Celery
 
Pre-requisites are: pip, virtualenv, PostgreSQL, RabbitMQ

From your virtual environment:

* Create a PostgreSQL database and update the database configuration in 
paperstream/config/settings/base.py with your username and database name

* Setup the project (dependencies, local files, migrations, initial data):  

```
#!bash

$ ./scripts/install_local.sh
```

* Fire-up the server:

```
#!bash
# Start RabbitMQ-server
$ rabbitmq-server

# Start Celery queues, and beat
$ celery -A config worker -Q default --concurrency=10 -l INFO -n worker1.default.%h
$ celery -A config worker -Q models,default --concurrency=10 -l INFO -n worker1.models.%h
$ celery -A config worker -Q consumers,default -l INFO -n worker1.consumers.%h
$ celery -A config beat -s ../logs/celerybeat-schedule

# Start WebApp
$ ./manage.py runserver

```

* Visit 127.0.0.1:8000

### Who do I talk to? ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)


