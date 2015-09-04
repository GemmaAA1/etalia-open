# -*- coding: utf-8 -*-
"""
Fabfile for paperstream deployment on AWS.

Current Configuration is as follow:
Two stacks are defined:
    - staging (at staging-stack.paperstream.io)
    - production (at www.paperstream.io)
Each stack has has 3 layers:
    - a Postgres db (host on AWS RDS)
    - apps
    - jobs

Each instance is tagged accordingly. For example, an instance with tags
{'site': 'staging', 'layer': 'jobs', 'Name': 'job1'} corresponds to an instance
that belongs to stack staging and layer jobs.

Deployment examples:
- deploy to layer jobs on staging stack:
    fab set_hosts:staging,jobs deploy
- deploy instance App1 on layer apps on the production stack:
    fab set_hosts:production,apps,App1 deploy
- deploy on all instances of stack staging:
    fab deploy:staging

NB: During first deployment, fabfile will stop to allow you to upload you id_rsa.pub
to Bitbucket to allow pulling

"""
from __future__ import absolute_import, unicode_literals, print_function
import os
from boto.ec2 import connect_to_region
from fabric.api import env, run, cd, settings, prefix, task, local, prompt, put, sudo
from fabric.contrib import files
from fabtools.utils import run_as_root
import fabtools

STACK = 'staging'
STACK_SITE_MAPPING = {'staging': 'staging-stack.paperstream.io',
                      'production': 'www.paperstream.io'}
REGION = os.environ.get("DJANGO_AWS_REGION")
SSH_EMAIL = 'nicolas.pannetier@gmail.com'
REPO_URL = 'git@bitbucket.org:NPann/paperstream.git'
VIRTUALENV_DIR = '.virtualenvs'
USER = 'ubuntu'


# Server user, normally AWS Ubuntu instances have default user "ubuntu"
# List of AWS private key Files
env.key_filename = ['~/.ssh/npannetier-key-pair-oregon.pem']
env.user = USER
env.virtualenv_dir = VIRTUALENV_DIR


def deploy(stack=STACK):
    # setup env
    setattr(env, 'stack', stack)
    setattr(env, 'stack_site', STACK_SITE_MAPPING.get(env.stack))
    setattr(env, 'env_dir', '/home/{user}/{venv}/{stack}'.format(
        user=env.user,
        venv=env.virtualenv_dir,
        stack=env.stack))
    if not env.hosts:
        set_hosts(stack=stack)
    stack_folder = '/home/{0}/{1}'.format(env.user, env.stack)
    setattr(env, 'site_folder', stack_folder)
    source_folder = '{0}/source'.format(stack_folder)
    setattr(env, 'source_folder', source_folder)
    # run
    _update_and_require_libraries()
    _create_virtual_env_if_necessary()
    _create_directory_structure_if_necessary()
    _get_latest_source()
    _pip_install()
    _update_static_files()
    _update_database()
    _update_nginx_conf_file()


def set_hosts(stack=STACK, layer='*', name='*', region=REGION):
    """Fabric task to set env.hosts based on tag key-value pair"""
    context = {
        'tag:stack': stack,
        'tag:layer': layer,
        'tag:Name': name,
    }
    env.hosts = _get_public_dns(region, context)


def _get_public_dns(region, context):
    """Private method to get public DNS name for instance with given tag key
     and value pair"""
    public_dns = []
    connection = _create_connection(region)
    reservations = connection.get_all_instances(filters=context)
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
                                   'rabbitmq-server',
                                   'supervisor',
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
        "workon {stack}".format(stack=env.stack),
        "export DJANGO_SETTINGS_MODULE=config.settings.{stack}".format(stack=env.stack) % env
    ]
    return prefix(" && ".join(workon_command))


def _create_directory_structure_if_necessary():
    for sub_folder in ('static', 'source'):
        if not files.exists('{0}/{1}'.format(env.stack_folder, sub_folder)):
            run('mkdir -p {0}/{1}'.format(env.stack_folder, sub_folder))


def _get_latest_source():
    # Generating public key for ssh bitbucket
    if not files.exists('/home/{}/.ssh/id_rsa.pub'.format(env.user)):
        print('Generate id_rsa for BitBucket git ssh\n')
        run_as_root('ssh-keygen')
        run_as_root('ps -e | grep [s]sh-agent')
        run_as_root('ssh-agent /bin/bash')
        run_as_root('ssh-add ~/.ssh/id_rsa ')
        print('Add the public below to your bitbutcket and run again:'
              '(https://confluence.atlassian.com/bitbucket/set-up-ssh-for-git-728138079.html)')
        run_as_root('cat ~/.ssh/id_rsa.pub')
        return
    if files.exists(env.source_folder + '/.git'):
        run('cd {0} && git fetch'.format(env.source_folder))
    else:
        run('git clone {0} {1}'.format(REPO_URL, env.source_folder))
    # git match what is checkout on your local repo
    current_commit = local("git log -n 1 --format=%H", capture=True)
    run('cd {0} && sudo git reset --hard {1}'.format(env.source_folder,
                                                     current_commit))


def _pip_install():
    with settings(cd(env.source_folder), _workon()):
        run('pip install -r requirements/{stack}.txt'.format(site=env.stack))


def _update_static_files():
    with settings(cd(env.source_folder), _workon()):
        run('cd {} && ./manage.py collectstatic --noinput'
            .format(env.source_folder))


def _update_database():
    with settings(cd(env.source_folder), _workon()):
        run('cd {} && ./manage.py migrate --noinput'.format(env.source_folder,))


def _update_nginx_conf_file():
    # remove file if exists
    available_file = '/etc/nginx/sites-available/{site}'.format(site=env.stack_site)
    if files.exists(available_file):
        run_as_root('rm ' + available_file)
    # upload template
    files.upload_template(
        '{source}/deploy/nginx.template.conf'.format(source=env.source_folder),
        available_file,
        context={'SITENAME': env.stack_site,
                 'USER': env.user,
                 'STACK': env.stack})
    # update link
    enable_link = '/etc/nginx/sites-enabled/{site}'.format(site=env.stack_site)
    if files.is_link(enable_link):
        run_as_root('rm {link}'.format(link=enable_link))
    run_as_root('ln -s /etc/nginx/sites-available/{site} '
                '/etc/nginx/sites-enabled/{site}'.format(site=env.stack_site))


def _update_gunicorn_start():
    # remove file if exists
    gunicorn_start_path = '{env_dir}/bin/gunicorn-start'.format(env_dir=env.env_dir)
    if files.exists(gunicorn_start_path):
        run_as_root('rm ' + gunicorn_start_path)
    # upload template
    files.upload_template(
        '{source}/deploy/gunicorn-start.template'.format(source=env.source_folder),
        gunicorn_start_path,
        context={'SITENAME': env.stack_site,
                 'SOURCE_FOLDER': env.source_folder,
                 'USER': env.user,
                 'STACK': env.stack})
    # change mod to execute
    run_as_root('sudo chmod 775 ' + gunicorn_start_path)