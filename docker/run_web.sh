#!/bin/sh

# prepare init migration
su -m myuser -c "python manage.py makemigrations etalia"
# migrate db, so we have the latest db schema
su -m myuser -c "python manage.py migrate"
# update etalia
su -m myuser -c "python setup/manager.py --update"
# start development server on public ip interface, on port 8000
su -m myuser -c "python manage.py runserver 0.0.0.0:8000"
