#!/usr/bin/env bash

# TODO Remove scripts from ./scripts/docker

# Create network if not exist
if [[ "" == "$(docker network ls | grep etalia-network)" ]]
then
    docker network create etalia-network
fi
# Create volume if not exist
if [[ "" == "$(docker volume ls | grep etalia-data-volume)" ]]
then
    docker volume create --name=etalia-data-volume
fi

Title() {
    printf "\n\e[1;104m ----- $1 ----- \e[0m\n"
}

Warning() {
    printf "\n\e[31;43m$1\e[0m\n"
}

Help() {
    printf "\n\e[2m$1\e[0m\n";
}

Confirm () {
    printf "\n"
    choice=""
    while [ "$choice" != "n" ] && [ "$choice" != "y" ]
    do
        printf "Do you want to continue ? (N/Y)"
        read choice
        choice=$(echo ${choice} | tr '[:upper:]' '[:lower:]')
    done
    if [ "$choice" = "n" ]; then
        printf "\nAbort by user.\n"
        exit 0
    fi
    printf "\n"
}

Compose() {
    docker-compose -f ./docker/docker-compose.yml $1
}

Run() {
    Compose "run worker $1"
}

Migrate() {
    Run "python manage.py makemigrations"
    Run "python manage.py migrate"
}

Load() {
    Run "python control/manager.py --load"
}

Update() {
    Run "python control/manager.py --update"
}

case $1 in
    build)
        docker build -f ./docker/Dockerfile -t etalia/dev .
    ;;
    install)
        Title "Installation"
        Confirm

        Compose "up -d"

        Migrate
        Load

        Compose "restart worker"

        Update
    ;;
    up)
        Compose "up -d"
    ;;
    down)
        Compose "down -v --remove-orphans"
    ;;
    dump)
        Run "python control/manager.py --dump"
    ;;
    init)
        Run "python control/manager.py --init"
    ;;
    load)
        Title "Load"
        Confirm

        Load
    ;;
    migrate)
        Title "Migrate"
        Confirm

        Migrate
    ;;
    update)
        Title "Update"
        Confirm

        Update
    ;;
    # ------------- HELP -------------
    *)
        printf "\n\e[2mUsage: ./setup.sh [action]

 - \e[0mbuild\e[2m\t Builds the container image.
 - \e[0minstall\e[2m\t Initial installation.
 - \e[0mup\e[2m\t\t Stack up.
 - \e[0mdown\e[2m\t\t Stack down.\e[0m\n"
    ;;
esac
