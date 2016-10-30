#!/bin/sh

# start development server on public ip interface, on port 8000
su -m myuser -c "python manage.py runserver 0.0.0.0:8000"
