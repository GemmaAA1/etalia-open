# SETUP DOCKER with ssh #

1. Build the docker image with ssh client:

        $ ./ssh/build-dev.sh

2. Run docker-compose for pycharm version

        $ docker-compose -f ./docker/docker-compose-pycharm.yml up -d

3. Configure a remote interpreter with ssh credentials in Pycharm:

    On MacOs:
        - Go to PyCharm/Preferences/Project/Project Interpreter
        - Add remote interpreter with settings as shown in attached screenshot

Note:
You can ssh to container based on python-dev-ssh running:

        ssh -i ./docker/ssh/ssh_keys/myuser-sshd -p 2222 myuser@localhost