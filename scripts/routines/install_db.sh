#!/usr/bin/env bash

# Cheat sheet for installing postgres

# Install postgres (with brew)
brew install postgres

# Init db with default username
initdb /usr/local/var/postgres -E utf8

# Creating a Database Cluster
initdb -D /usr/local/pgsql/data

# Starting database
pg_ctl -D /usr/local/pgsql/data -l logfile start

# create <user> database
createdb

