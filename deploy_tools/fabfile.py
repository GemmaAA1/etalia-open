# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals, print_function
import os
from boto.ec2 import connect_to_region
from fabric.api import env, run, cd, settings, prefix, task, local, prompt, put, sudo
from fabric.contrib import files
from fabtools.utils import run_as_root
import fabtools

SITE = 'staging'
REGION = os.environ.get("DJANGO_AWS_REGION")
SSH_EMAIL = 'nicolas.pannetier@gmail.com'
REPO_URL = 'git@bitbucket.org:NPann/paperstream.git'
VIRTUALENV_DIR = '.virtualenvs'

# Server user, normally AWS Ubuntu instances have default user "ubuntu"
env.user = 'ubuntu'
env.site = SITE
# List of AWS private key Files
env.key_filename = ['~/.ssh/npannetier-key-pair-oregon.pem']
env.virtualenv_dir = VIRTUALENV_DIR


def deploy(site=SITE):
    # setup
    env.site = site
    if not env.hosts:
        set_hosts(tag=site, value='*')
    site_folder = '/home/{0}/{1}'.format(env.user, env.site)
    setattr(env, 'site_folder', site_folder)
    source_folder = '{0}/source'.format(site_folder)
    setattr(env, 'source_folder', source_folder)
    # run
    _update_and_require_libraries()
    _create_virtual_env_if_necessary()
    _create_directory_structure_if_necessary()
    _get_latest_source()
    _pip_install()
    _update_static_files()
    # _update_database(source_folder)


def set_hosts(tag=SITE, value='*', region=REGION):
    """Fabric task to set env.hosts based on tag key-value pair"""
    key = "tag:"+tag
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
    fabtools.require.deb.packages(['python3-dev',
                                   'build-essential',
                                   'python3-setuptools',
                                   'python3-numpy',
                                   'python3-scipy',
                                   'python3-pip',
                                   'libblas-dev',
                                   'liblapack-dev',
                                   'libatlas-dev',
                                   'libatlas3gf-base',
                                   'gfortran',
                                   'libncurses5-dev',
                                   'postgresql',
                                   'nginx',
                                   'git',
                                   'python-pip',
                                   'libpq-dev',
                                   ], update=True)
    # Require some pip packages
    fabtools.require.python.packages(['virtualenvwrapper', ], use_sudo=True)
    # Set virtualenvwrapper
    if not files.contains("~/.bashrc", "export WORKON_HOME={virtualenv_dir}".format(
            virtualenv_dir=env.virtualenv_dir)):
        files.append("~/.bashrc", "export WORKON_HOME={virtualenv_dir}".format(
            virtualenv_dir=env.virtualenv_dir))
    if not files.contains("~/.bashrc", "source /usr/local/bin/virtualenvwrapper.sh"):
        files.append("~/.bashrc", "source /usr/local/bin/virtualenvwrapper.sh")
    if not files.exists(env.virtualenv_dir):
        run("mkdir -p {virtualenv_dir}".format(virtualenv_dir=env.virtualenv_dir))


def _create_virtual_env_if_necessary():
    # Create virtual env
    with prefix("WORKON_HOME={virtualenv_dir}".format(virtualenv_dir=env.virtualenv_dir)):
        with prefix('source /usr/local/bin/virtualenvwrapper.sh'):
            existent_virtual_envs = run("lsvirtualenv")
            if not env.site in existent_virtual_envs:
                run("mkvirtualenv --python=/usr/bin/python3 {virtual_env}".format(virtual_env=env.site))


def _workon():
    workon_command = [
        "source /usr/local/bin/virtualenvwrapper.sh",
        "workon {site}".format(site=env.site),
        "export DJANGO_SETTINGS_MODULE=config.settings.{site}".format(site=env.site) % env
    ]
    return prefix(" && ".join(workon_command))


def _create_directory_structure_if_necessary():
    for sub_folder in ('static', 'source'):
        if not files.exists('{0}/{1}'.format(env.site_folder, sub_folder)):
            run('mkdir -p {0}/{1}'.format(env.site_folder, sub_folder))


def _get_latest_source(source_folder):
    # Generating public key for ssh bitbucket
    if not files.exists('/home/{}/.ssh/id_rsa.pub'.format(env.user)):
        print('Generate id_rsa for BitBucket git ssh\n')
        run_as_root('ssh-keygen')
        run_as_root('ps -e | grep [s]sh-agent')
        run_as_root('ssh-agent /bin/bash')
        run_as_root('ssh-add ~/.ssh/id_rsa ')
        print('Add the public below to your bitbutcket and run again'
              '(https://confluence.atlassian.com/bitbucket/set-up-ssh-for-git-728138079.html)')
        run('cat ~/.ssh/id_rsa.pub')
        return
    if files.exists(source_folder + '/.git'):
        run('cd {0} && git fetch'.format(env.source_folder))
    else:
        run('git clone {0} {1}'.format(REPO_URL, env.source_folder))
    # git match what is checkout on your local repo
    current_commit = local("git log -n 1 --format=%H", capture=True)
    run('cd {0} && sudo git reset --hard {1}'.format(env.source_folder,
                                                     current_commit))


def _pip_install():
    with settings(cd(env.source_folder), _workon()):
        run('pip install -r requirements/{site}.txt'.format(site=env.site))


def _update_static_files():
    with settings(cd(env.source_folder), _workon()):
        run('cd {} && ./manage.py collectstatic --noinput'
            .format(env.source_folder))

