#!/bin/sh

docker build -f ./docker/Dockerfile -t etalia/python-dev .

# Generate or re-use SSH keys. The ssh_keys directory should be in .gitignore
mkdir -p ssh_keys
if [ ! -f ssh_keys/docker-sshd ] ; then
    echo "Generating SSH key"
    # Remove any associated public key also
    rm -f ssh_keys/docker-sshd*
    ssh-keygen -t rsa -f ssh_keys/docker-sshd -N ""
else
    echo "Using existing SSH key"
fi

# unfortunately docker native doesn't mount in /Applications, either change that or do this
rm -rf ~/.pycharm_helpers
cp -R /Applications/PyCharm.app/Contents/helpers ~/.pycharm_helpers

# Now build the sshd server image, which should be based on the main image.
docker build --force-rm=true -f docker/Dockerfile.ssh \
 --tag etalia/python-dev-ssh .

#--build-arg AUTH_KEY="$(cat ssh_keys/docker-sshd.pub)" \