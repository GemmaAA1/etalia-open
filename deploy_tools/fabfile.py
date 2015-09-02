# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals, print_function
import os
from boto.ec2 import connect_to_region
from fabric.api import *
from fabtools.python import virtualenv
from fabtools import require
from fabtools import files
from fabtools.utils import run_as_root, env
import fabtools

SITE = 'staging'

REGION = os.environ.get("DJANGO_AWS_REGION")
SSH_EMAIL = 'nicolas.pannetier@gmail.com'
REPO_URL = 'git@bitbucket.org:NPann/paperstream.git'
# Server user, normally AWS Ubuntu instances have default user "ubuntu"
env.user = "ubuntu"
# List of AWS private key Files
env.key_filename = ["~/.ssh/npannetier-key-pair-oregon.pem"]


def deploy(site=SITE):
    setattr(env, 'site', site)
    if not env.hosts:
        set_hosts(tag=site)
    site_folder = '/home/{0}/{1}'.format(env.user, env.site)
    source_folder = '{0}/source'.format(site_folder)
    _update_and_require_libraries()
    _create_directory_structure_if_necessary(site_folder)
    _get_latest_source(source_folder)
    _update_virtualenv(source_folder)
    # _update_static_files(source_folder)
    # _update_database(source_folder)


def set_hosts(tag=SITE, value="*", region=REGION):
    """Fabric task to set env.hosts based on tag key-value pair"""
    key = "tag:"+tag
    env.site_folder = SITE
    env.hosts = _get_public_dns(region, key, value)


def _get_public_dns(region, key, value="*"):
    """Private method to get public DNS name for instance with given tag key
     and value pair"""
    public_dns = []
    connection = _create_connection(region)
    reservations = connection.get_all_instances(filters={key: value})
    for reservation in reservations:
        for instance in reservation.instances:
            print("Instance {}".format(instance.public_dns_name))
            public_dns.append(str(instance.public_dns_name))
    return public_dns


def _create_connection(region):
    """Private method for getting AWS connection"""
    print("Connecting to {}".format(region))
    conn = connect_to_region(
        region_name=region,
        aws_access_key_id=os.environ.get("DJANGO_AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.environ.get("DJANGO_AWS_SECRET_ACCESS_KEY")
    )
    print("Connection with AWS established")
    return conn


def _update_and_require_libraries():
    # Require some Ubuntu packages
    require.deb.packages(['python-dev',
                          'postgresql',
                          'postgresql-contrib',
                          'nginx',
                          'git',
                          'python-pip',
                          'python-psycopg2',
                          'libpq-dev',
                          'libncurses5-dev',
                          'libatlas-base-dev',
                          'python-numpy',
                          'python-scipy',
                          'gfortran',
                          'rabbitmq-server',
                          'libpq-dev',
                          ], update=True)
    # Require some pip packages
    require.python.packages(['virtualenv', ],  use_sudo=True)


def _create_directory_structure_if_necessary(site_folder):
    for sub_folder in ('static', 'source'):
        if not files.is_dir('{0}/{1}'.format(site_folder, sub_folder)):
            run('mkdir -p {0}/{1}'.format(site_folder, sub_folder))


def _get_latest_source(source_folder):
    # Generating public key for ssh bitbucket
    if not files.is_file('/home/{}/.ssh/id_rsa.pub'.format(env.user)):
        print('Generate id_rsa for BitBucket git ssh\n')
        run('ssh-keygen')
        run('ps -e | grep [s]sh-agent')
        run('ssh-agent /bin/bash')
        run('ssh-add ~/.ssh/id_rsa ')
        print('Add the public below to your bitbutcket and run again'
              '(https://confluence.atlassian.com/bitbucket/set-up-ssh-for-git-728138079.html)')
        run('cat ~/.ssh/id_rsa.pub')
        return
    if files.is_dir(source_folder + '/.git'):
        run('cd {0} && git fetch'.format(source_folder))
    else:
        run('git clone {0} {1}'.format(REPO_URL, source_folder))
    # git match what is currently checked out on your local repo
    current_commit = local("git log -n 1 --format=%H", capture=True)
    run('cd {0} && sudo git reset --hard {1}'.format(source_folder,
                                                     current_commit))


def _update_virtualenv(source_folder):
    virtualenv_folder = source_folder + '/../virtualenv'
    if not files.is_file(virtualenv_folder + '/bin/pip'):
        run('virtualenv --python=python2.7 {}'.format(virtualenv_folder))
    # install requirement
    run('{0}/bin/pip install --upgrade pip'.format(virtualenv_folder))
    run('{0}/bin/pip install -r {1}/requirements/{2}.txt'.format(
        virtualenv_folder, source_folder, env.site))


