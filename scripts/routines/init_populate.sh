#!/usr/bin/env bash

# Check if manage.py exists
if [ ! -f '../manage.py' ]; then
    echo '../manage.py not found'
fi

## Populate database
# populate library
../manage.py populate publisher all
../manage.py populate journal thomson
../manage.py populate journal pubmed
../manage.py populate journal arxiv

# populate consumers
../manage.py populate consumer pubmed --name pubmed_all
../manage.py populate consumer arxiv --name arxiv_all
../manage.py populate consumer elsevier --name elsevier_all