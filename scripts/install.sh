#!/usr/bin/env bash

# Set environment variable for current VIRTUAL_ENV
./routines/init_env_vars.sh

# Create Postgres database
./routines/create_db.sh

# Init migrations and migrate
./routines/init_migrations.sh

# Populate database
./routines/init_populate.sh

