#!/usr/bin/env bash

# Check if virtual env is active
./routines/init_env_vars.sh

# Init migrations and migrate
./routines/init_migrations.sh

## Populate database
./routines/init_populate.sh

