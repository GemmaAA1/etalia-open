# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

The steps below will get you up and running with a local development environment. We assume you have the following installed:

pip
virtualenv
PostgreSQL
First make sure to create and activate a virtualenv, then open a terminal at the project root and install the requirements for local development:
```
#!bash

$ pip install -r requirements/local.txt
```

Then, create a PostgreSQL database and add the database configuration in the config.settings.base.py setting file

You can now run the usual Django migrate and runserver command:
```
#!bash

$ python manage.py migrate
```
```
#!bash

$ python manage.py runserver
```
### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact