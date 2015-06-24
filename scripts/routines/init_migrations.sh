#!/usr/bin/env bash

# Check if manage.py exists
if [ ! -f '../paperstream/manage.py' ]; then
    echo '../paperstream/manage.py not found'
fi

../paperstream/manage.py makemigrations
../paperstream/manage.py migrate