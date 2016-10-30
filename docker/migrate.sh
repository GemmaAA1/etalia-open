#!/bin/sh
# prepare init migration
python manage.py makemigrations
# migrate db, so we have the latest db schema
python manage.py migrate