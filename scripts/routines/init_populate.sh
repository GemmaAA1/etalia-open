#!/usr/bin/env bash

# Check if manage.py exists
if [ ! -f '../paperstream/manage.py' ]; then
    echo '../paperstream/manage.py not found'
fi

## Populate database
# populate library
../paperstream/manage.py populate publisher all
../paperstream/manage.py populate journal thomson
../paperstream/manage.py populate journal pubmed
../paperstream/manage.py populate journal arxiv

# populate consumers
../paperstream/manage.py populate consumer pubmed --name pubmed_all
../paperstream/manage.py populate consumer arxiv --name arxiv_all
../paperstream/manage.py populate consumer elsevier --name elsevier_all