#!/bin/sh

CUR_DIR="$(pwd)"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname ""$(dirname "$DIR")"")"
cd $DIR

# Generate or re-use SSH keys. The ssh_keys directory should be in .gitignore
mkdir -p ssh_keys
if [ ! -f ssh_keys/myuser-sshd ] ; then
    echo "Generating SSH key"
    # Remove any associated public key also
    rm -f ssh_keys/myuser-sshd*
    ssh-keygen -t rsa -f ssh_keys/myuser-sshd -N ""
else
    echo "Using existing SSH key"
fi

# NB: This is only to use remote interpreter from pycharm within the docker container
# unfortunately docker native doesn't mount in /Applications, either change that or do this
rm -rf ~/.pycharm_helpers
cp -R /Applications/PyCharm.app/Contents/helpers ~/.pycharm_helpers

# Now build the sshd server image, which should be based on the main image.
cd $ROOT_DIR
docker build --force-rm=true -f docker/ssh/Dockerfile --tag etalia/dev-ssh .

cd $CUR_DIR