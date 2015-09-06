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
    >> fab set_hosts:staging,jobs deploy
- deploy instance app1 on layer apps on the production stack:
    >> fab set_hosts:production,apps,app1 deploy
- reset env var
    >> fab set_hosts:staging,apps,app1 set_virtual_env_hooks

NB:
- During first deployment, fabfile will stop to allow you to upload you id_rsa.pub
to Bitbucket to allow pulling
- Environment variables are defined

"""
from __future__ import absolute_import, unicode_literals, print_function
import os
import re
import string
import random
from boto.ec2 import connect_to_region
from fabric.decorators import roles, runs_once, task, parallel
from fabric.api import env, run, cd, settings, prefix, task, local, prompt, put, sudo, get
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
SUPERVISOR_CONF_DIR = 'supervisor'
USER = 'ubuntu'


# Server user, normally AWS Ubuntu instances have default user "ubuntu"
# List of AWS private key Files
env.key_filename = ['~/.ssh/npannetier-key-pair-oregon.pem']
env.user = USER
env.virtualenv_dir = VIRTUALENV_DIR
env.conf_dir = SUPERVISOR_CONF_DIR

@task
def set_hosts(stack=STACK, layer='*', name='*', region=REGION):
    """Fabric task to set env.hosts based on tag key-value pair"""
    # setup env
    setattr(env, 'stack', stack)
    setattr(env, 'stack_site', STACK_SITE_MAPPING.get(env.stack))
    setattr(env, 'env_dir', '/home/{user}/{venv}/{stack}'.format(
        user=env.user,
        venv=env.virtualenv_dir,
        stack=env.stack))
    stack_dir = '/home/{0}/{1}'.format(env.user, env.stack)
    setattr(env, 'stack_dir', stack_dir)
    source_dir = '{0}/source'.format(stack_dir)
    setattr(env, 'source_dir', source_dir)

    context = {
        'tag:stack': stack,
        'tag:layer': layer,
        'tag:Name': name,
    }
    env.hosts, env.roles, env.roledefs = _get_public_dns(region, context)


def get_host_roles():
    return [k for k, v in env.roledefs.items() if env.host_string in v]


@task
def deploy():
    """Deploy paperstream on hosts"""
    if not env.hosts:
        raise ValueError('No hosts defined')
    # run
    update_and_require_libraries()
    create_virtual_env_if_necessary()
    create_virtual_env_hooks_if_necessary()
    create_directory_structure_if_necessary()
    pull_latest_source()
    pip_install()
    update_database()
    update_static_files()
    if env.host_string in env.roledefs.get('apps', []):
        set_rabbit_user()
        update_gunicorn_conf()
    if env.host_string in env.roledefs.get('jobs', []):
        update_supervisor_conf()
        update_nginx_conf()
        reload_nginx()
    reload_supervisor()


def create_virtual_env_hooks_if_necessary():
    """Create environment variable hooks if not defined"""
    postactivate_path = '{env_dir}/bin/postactivate'.format(env_dir=env.env_dir)
    if not files.contains(postactivate_path, 'DJANGO_SECRET_KEY'):  # this is the bare minimum
        update_virtual_env_hooks()


@task
def update_virtual_env_hooks():
    """Define environment variable on host based on local file"""
    if not env.hosts:
        raise ValueError('No hosts defined')

    # remove file if exists
    postactivate_path = '{env_dir}/bin/postactivate'.format(env_dir=env.env_dir)
    if files.exists(postactivate_path):
        run_as_root('rm ' + postactivate_path)
    # update file base on local file
    run('touch ' + postactivate_path)
    run('chmod 644 ' + postactivate_path)
    run('echo "#!/usr/bin/env bash" >>' + postactivate_path)
    env_vars = []
    with open('.envs/{stack}'.format(stack=env.stack)) as f:  # on local machine
        for line in f:
            if line.strip():  # prevent from empty line
                if 'DJANGO_SECRET_KEY' in line:  # generate key
                    key = _generate_new_secret_key(key_length=32)
                    line = 'DJANGO_SECRET_KEY="{key}"'.format(key=key)
                run("echo 'export {line}' >> {dest}".format(
                    line=line.strip(), dest=postactivate_path))
                # store env var
                res = re.match(r'([A-Z0-9_]+)=', line)
                if res.groups():
                    env_vars.append(res.groups()[0])

    # update predeactivate
    predeactivate_path = '{env_dir}/bin/predeactivate'.format(env_dir=env.env_dir)
    if files.exists(predeactivate_path):
        run_as_root('rm ' + predeactivate_path)
    run('touch ' + predeactivate_path)
    run('chmod 644 ' + predeactivate_path)
    run('echo "#!/usr/bin/env bash" >>' + predeactivate_path)
    for env_var in env_vars:
        run("echo 'unset {var}' >> {dest}".format(
            var=env_var, dest=predeactivate_path))


def _get_public_dns(region, context):
    """Private method to get public DNS name for instance with given tag key
     and value pair"""
    public_dns = []
    tags = {}
    connection = _create_connection(region)
    reservations = connection.get_all_instances(filters=context)
    for reservation in reservations:
        for instance in reservation.instances:
            if instance.public_dns_name:
                print("Instance {}".format(instance.public_dns_name))
                public_dns.append(str(instance.public_dns_name))
                tags[str(instance.public_dns_name)] = instance.tags

    # Define rol    es for instances based on tag
    roles = list(set([tags[host]['layer'] for host in public_dns]))
    roledefs = {}
    for role in roles:
        roledefs[role] = [host for host in public_dns
                          if tags[host]['layer'] == role]

    return public_dns, roles, roledefs


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


def _generate_new_secret_key(key_length=32):
    """Return new secret key"""
    symbols = string.ascii_letters + str(string.digits) + '!@#$%^&*()+-'
    secret_key = "".join(random.sample(symbols * 2, key_length))
    return secret_key


@task
def update_and_require_libraries():
    """Update ubuntu libraries"""
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
                                   'supervisor',
                                   ], update=True)
    # Require some pip packages
    fabtools.require.python.packages(['virtualenvwrapper'], use_sudo=True)

    # Set virtualenvwrapper
    if not files.contains("~/.bashrc", "export WORKON_HOME={virtualenv_dir}".format(
            virtualenv_dir=env.virtualenv_dir)):
        files.append("~/.bashrc", "export WORKON_HOME={virtualenv_dir}".format(
            virtualenv_dir=env.virtualenv_dir))
    if not files.contains("~/.bashrc", "source /usr/local/bin/virtualenvwrapper.sh"):
        files.append("~/.bashrc", "source /usr/local/bin/virtualenvwrapper.sh")
    if not files.exists(env.virtualenv_dir):
        run("mkdir -p {virtualenv_dir}".format(virtualenv_dir=env.virtualenv_dir))



def create_virtual_env_if_necessary():
    # Create virtual env
    with prefix("WORKON_HOME={virtualenv_dir}".format(virtualenv_dir=env.virtualenv_dir)):
        with prefix('source /usr/local/bin/virtualenvwrapper.sh'):
            existent_virtual_envs = run("lsvirtualenv")
            if not env.stack in existent_virtual_envs:
                run("mkvirtualenv --python=/usr/bin/python3 {virtual_env}".format(virtual_env=env.stack))


def _workon():
    workon_command = [
        "source /usr/local/bin/virtualenvwrapper.sh",
        "workon {stack}".format(stack=env.stack),
        "export DJANGO_SETTINGS_MODULE=config.settings.{stack}".format(stack=env.stack) % env
    ]
    return prefix(" && ".join(workon_command))


@task
def create_directory_structure_if_necessary():
    for sub_dir in ('static', 'source', 'log', env.conf_dir):
        if not files.exists('{0}/{1}'.format(env.stack_dir, sub_dir)):
            run('mkdir -p {0}/{1}'.format(env.stack_dir, sub_dir))


@task
def pull_latest_source():
    """Pull source from bitbucket"""
    # Test if private key has been uploaded
    if not files.exists('/home/{}/.ssh/id_rsa'.format(env.user)):
        if not files.exists('/home/{}/.ssh'.format(env.user)):
            run('mkdir /home/{}/.ssh'.format(env.user))
        put('.ssh/id_rsa', '/home/{}/.ssh/id_rsa'.format(env.user))
        run('chmod 400 /home/{}/.ssh/id_rsa'.format(env.user))
    # pull git
    if files.exists(env.source_dir + '/.git'):
        run('cd {0} && git fetch'.format(env.source_dir))
    else:
        run('git clone {0} {1}'.format(REPO_URL, env.source_dir))
    # git match what is checkout on your local repo
    current_commit = local("git log -n 1 --format=%H", capture=True)
    run('cd {0} && sudo git reset --hard {1}'.format(env.source_dir,
                                                     current_commit))

@task
def pip_install():
    """Pip install requirements"""
    with settings(cd(env.source_dir), _workon()):
        run('pip install -r requirements/{stack}.txt'.format(stack=env.stack))

@task
def update_static_files():
    """Update static files"""
    with settings(cd(env.source_dir), _workon()):
        run('cd {} && ./manage.py collectstatic --noinput'
            .format(env.source_dir))

@task
@runs_once
def update_database():
    """Update database from migrations"""
    with settings(cd(env.source_dir), _workon()):
        run('cd {} && ./manage.py migrate --noinput'.format(env.source_dir,))


@task
@roles('apps')
def update_nginx_conf():
    """Update Nginx conf files"""
    # remove file if exists
    available_file = '/etc/nginx/sites-available/{site}'.format(site=env.stack_site)
    if files.exists(available_file):
        run_as_root('rm ' + available_file)
    # upload template
    files.upload_template('nginx.template.conf', available_file,
        context={'SITENAME': env.stack_site,
                 'USER': env.user,
                 'STACK': env.stack}, use_sudo=True, use_jinja=True)
    # update link
    enable_link = '/etc/nginx/sites-enabled/{site}'.format(site=env.stack_site)
    if files.is_link(enable_link):
        run_as_root('rm {link}'.format(link=enable_link))
    run_as_root('ln -s /etc/nginx/sites-available/{site} '
                '/etc/nginx/sites-enabled/{site}'.format(site=env.stack_site))

@task
@roles('apps')
def update_gunicorn_conf():
    """Update Gunicorn conf files"""
    # remove file if exists
    gunicorn_conf_path = '{stack_dir}/{conf_dir}/gunicorn.conf.py'.format(
        stack_dir=env.stack_dir,
        conf_dir=env.conf_dir)
    if files.exists(gunicorn_conf_path):
        run_as_root('rm ' + gunicorn_conf_path)
    # upload template
    files.upload_template( 'gunicorn.template.conf.py', gunicorn_conf_path,
        context={'SITENAME': env.stack_site,
                 'USER': env.user}, use_sudo=True, use_jinja=True)


@task
@roles('jobs')
def set_rabbit_user():
    """Set rabbit user"""
    with settings(_workon()):
        if not files.exists("/usr/sbin/rabbitmq-server"):
            fabtools.require.deb.packages(['rabbitmq-server'])
            run_as_root('rabbitmqctl add_user $RABBITMQ_USERNAME $RABBITMQ_PASSWORD')
            run_as_root('rabbitmqctl set_permissions $RABBITMQ_USERNAME ".*" ".*" ".*"')
            run_as_root("rabbitmqctl delete_user guest")


@task
def update_supervisor_conf():
    """Set supervisor conf file"""
    if not files.exists('/var/log/supervisord'):
        run_as_root('mkdir -p /var/log/supervisord')

    # remove file if exists
    supervisor_file = '{stack_dir}/{conf_dir}/supervisord.conf'.format(
        stack_dir=env.stack_dir,
        conf_dir=env.conf_dir)
    if files.exists(supervisor_file):
        run_as_root('rm ' + supervisor_file)
    # upload template
    files.upload_template('supervisord.template.conf', supervisor_file,
        context={'SITENAME': env.stack_site,
                 'USER': env.user,
                 'STACK': env.stack,
                 'STACK_DIR': env.stack_dir,
                 'SOURCE_DIR': env.source_dir,
                 'ENV_DIR': env.env_dir,
                 'CONF_DIR': env.conf_dir,
                 'ROLES': get_host_roles(),
                 }, use_sudo=True, use_jinja=True)

    # Copy env variable from postactivate
    run('python {source_dir}/deploy/cp_p2s.py -i {env_dir}/bin/postactivate -o {supervisor}'.format(
        source_dir=env.source_dir,
        env_dir=env.env_dir,
        supervisor=supervisor_file))

    # simlink
    if files.exists('/etc/supervisor/conf.d/supervisord.conf'):
        run_as_root('rm /etc/supervisor/conf.d/supervisord.conf')
    run_as_root('ln -s {stack_dir}/{conf_dir}/supervisord.conf /etc/supervisor/conf.d/supervisord.conf'.format(
        stack_dir=env.stack_dir,
        conf_dir=env.conf_dir))


@task
def reload_nginx():
    run_as_root('sudo service nginx reload')


@task
def reload_supervisor():
    run_as_root('supervisorctl reread')
    run_as_root('supervisorctl update')

