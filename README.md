# README #

PaperStream is a WebApp for researchers and others who are reading scientific articles.

### What is PaperStream ? ###

* A personalized stream of scientific articles 
* A commenting tool to bring papers alive
* A way to explore the scientific literature

### How do I get set up? ###

PaperStream is developed with Django/PostgreSQL. Asynch tasks are managed with
RabbitMQ and Celery.  

The steps below will get you up and running with a local development environment. 
Pre-requisites installs are: pip, virtualenv, PostgreSQL, RabbitMQ

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
# Start Celery
$ celery -A config worker --loglevel=info
# Start WebApp
$ ./manage.py runserver

```

* Visit 127.0.0.1:8000

### Who do I talk to? ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)
* Valentine Toulemonde


