# README #

Etalia is a web application that facilitates scientific communications through 
text-driven technology 

### How do I get set up locally? ###

Etalia stack is based on: 

* Django 
* PostgreSQL
* RabbitMQ
* Celery

Requirements are: 

* pip (the python package manager)
* virtualenv
* PostgreSQL
* RabbitMQ

Once you get the requirements installed, follow these steps to get started:

1. Create you virtual environment and export environment variables from the 
development configuration file to your virtual environment. 
(Request one if you don't have one already)
2. Clone Etalia repository
```
#!bash

git clone https://NPann@bitbucket.org/NPann/etalia.git
```

3. Create a local Postgres database (e.g. "CREATE DATABASE etalia;" from psql on mac)
4. Copy etalia/config/settings/common.py.dist to etalia/config/settings/common.py
4. Update database settings in etalia/config/settings/common.py with your 
username and database name. Example:
 
```
DATABASES = {
  'default': {
      ...
      'NAME': 'database_name',
      'USER': 'your_username',
      'PASSWORD': '',
      ...
  }
}
```
   
3. Go to scripts directory and run the install_local script:
```
cd scripts/
./install_local.sh
``` 


* Launching the app:

```
$ ./manage.py runserver
```

* Visit 127.0.0.1:8000

### Etalia team ###

* Nicolas Pannetier (nicolas.pannetier@gmail.com)
* Etienne Dauvergne
* Norbert Schuff

