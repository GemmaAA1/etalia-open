# README #

PaperStream is a WebApp for researchers and others who are reading scientific articles.

### What is PaperStream ? ###

* A stream of curated articles 
* A tool to bring paper alive
* A way to explore the scientific literature
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

PaperStream is developed with Django/PostgreSQL. Asynch tasks are managed with
RabbitMQ and Celery.  

The steps below will get you up and running with a local development environment. 
Pre-requisites installs are:
* pip
* virtualenv
* PostgreSQL
* RabbitMQ

Once your virtual environment is created, activate it and:

1. Create a PostgreSQL database and update the database configuration in 
paperstream/config/settings/base.py with your username and database name

2. Setup the project (dependencies, local files, migrations, initial data):  

```
#!bash

$ ./scripts/install_local.sh
```

3. Fire-up the server:

```
#!bash

$ ./manage.py runserver
# Start Celery
$ celery -A config worker --loglevel=info
# Start RabbitMQ-server
$ rabbitmq-server

```

4. Visit 127.0.0.1:8000

### Who do I talk to? ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)
* Valentine Toulemonde


