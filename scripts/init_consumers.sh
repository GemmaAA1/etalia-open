#!/bin/bash

./manage.py populate consumer pubmed --name pubmed_all
./manage.py populate consumer arxiv --name arxiv_all
./manage.py populate consumer elsevier --name elsevier_all