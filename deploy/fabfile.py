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
    - jobs (one master (JOB_MASTER) and satellites)

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
from fabric.api import env, run, cd, settings, prefix, task, local, prompt, \
    put, sudo, get, reboot
from fabric.contrib import files
from fabtools.utils import run_as_root
import fabtools

# STACK = 'staging'
STACK_SITE_MAPPING = {'staging': 'alpha-u6VcayvtcI.pubstream.io',
                      'production': 'www.pubstream.io'}
REGION = os.environ.get("DJANGO_AWS_REGION")
SSH_EMAIL = 'nicolas.pannetier@gmail.com'
REPO_URL = 'git@bitbucket.org:NPann/paperstream.git'
VIRTUALENV_DIR = '.virtualenvs'
SUPERVISOR_CONF_DIR = 'supervisor'
USER = 'ubuntu'
JOB_MASTER = 'job1'
INSTANCE_TYPES_RANK = { 't2.micro': 0,
                        't2.small': 1,
                        't2.medium': 2,
                        't2.large': 3,
                        'm4.large': 3,
                        'm4.xlarge': 4,
                        'm4.2xlarge': 5,
                        'm4.10xlarge': 6}
ROLE_INSTANCE_TYPE_MAP = {'web': 't2.micro',
                          'base': 't2.small',
                          'ms': 't2.medium',
                          'feed': 't2.medium',
                          'nlp': 't2.large',
                          'spot': 'm4.large'}

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
    tags = _get_public_dns(region, context)

    env.hosts = tags.keys()
    roles = []
    roles += [tag.get('layer', '') for tag in tags.values()]
    for tag in tags.values():
        if tag.get('layer', None):
            roles.append(tag.get('layer'))
        if tag.get('role', None):
            r = tag.get('role', '').split(',')
            for rr in r:
                roles.append(rr)
    env.roles = list(set(roles))
    roledefs = {}
    for role in env.roles:
        roledefs[role] = [host for host in env.hosts
                          if tags[host]['layer'] == role or
                          role in tags[host].get('role', '')]
    env.roledefs = roledefs

    # store stack used
    setattr(env, 'stack_string', stack)
    # store tags
    setattr(env, 'tags', tags)

    check_integrity()


def check_integrity():
    """Check if instance type is compatible with role"""
    for host in env.hosts:
        inst_type = env.tags[host]['type']
        roles = env.tags[host].get('role', '').split(',')
        if roles:
            for role in roles:
                min_type = ROLE_INSTANCE_TYPE_MAP[role]
                if not INSTANCE_TYPES_RANK[min_type] <= INSTANCE_TYPES_RANK[inst_type]:
                    raise TypeError('Instance {0} too small, type is {1} (min {2})'.format(
                        host,
                        inst_type,
                        min_type
                    ))


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
    # app related
    if env.host_string in env.roledefs.get('apps', []):
        update_gunicorn_conf()
        update_nginx_conf()
        reload_nginx()
    # job related
    if env.host_string in env.roledefs.get('jobs', []):
        update_rabbit_user()
    update_supervisor_conf()
    restart_supervisor()
    reb = update_hosts_file(env.stack_string)
    if reb:
        reboot_instance()


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
                tags[str(instance.public_dns_name)]['type'] = instance.instance_type

    return tags


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
                                   'libjpeg-dev',
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


@task
def create_virtual_env_if_necessary():
    # Create virtual env
    with prefix("WORKON_HOME={virtualenv_dir}".format(virtualenv_dir=env.virtualenv_dir)):
        with prefix('source /usr/local/bin/virtualenvwrapper.sh'):
            existent_virtual_envs = run("lsvirtualenv")
            if not env.stack in existent_virtual_envs:
                run("mkvirtualenv --python=/usr/bin/python3 --system-site-packages {virtual_env}".format(virtual_env=env.stack))


@task
def remove_virtual_env():
    with prefix('source /usr/local/bin/virtualenvwrapper.sh'):
        run("rmvirtualenv {virtual_env}".format(virtual_env=env.stack))


def _workon():
    workon_command = [
        "source /usr/local/bin/virtualenvwrapper.sh",
        "workon {stack}".format(stack=env.stack),
        "export DJANGO_SETTINGS_MODULE=config.settings.{stack}".format(stack=env.stack) % env
    ]
    return prefix(" && ".join(workon_command))


@task
def create_directory_structure_if_necessary():
    for sub_dir in ('static', 'source', 'log', 'db', env.conf_dir):
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
def update_rabbit_user():
    """Set rabbit user on job_master"""
    if env.tags[env.host_string].get('Name') == JOB_MASTER:
        with settings(_workon()):  # to get env var
            # if rabbitmq not installed, install
            if not files.exists("/usr/sbin/rabbitmq-server"):
                fabtools.require.deb.packages(['rabbitmq-server'])
            # test if user exists:
            list_users = run_as_root('rabbitmqctl list_users')
            user = run_as_root('echo $RABBITMQ_USERNAME')
            if user not in list_users:
                run_as_root('rabbitmqctl add_user $RABBITMQ_USERNAME $RABBITMQ_PASSWORD')
                run_as_root('rabbitmqctl set_user_tags $RABBITMQ_USERNAME administrator')
                run_as_root('rabbitmqctl set_permissions $RABBITMQ_USERNAME ".*" ".*" ".*"')
            if 'guest' in list_users:
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
    with settings(_workon()):  # to get env var
        flower_users_passwords = run('echo $USERS_PASSWORDS_FLOWER')
        if env.tags[env.host_string].get('Name') == JOB_MASTER:
            job_master = True
        else:
            job_master = False
        files.upload_template('supervisord.template.conf', supervisor_file,
            context={'SITENAME': env.stack_site,
                     'USER': env.user,
                     'STACK': env.stack,
                     'STACK_DIR': env.stack_dir,
                     'SOURCE_DIR': env.source_dir,
                     'ENV_DIR': env.env_dir,
                     'CONF_DIR': env.conf_dir,
                     'ROLES': get_host_roles(),
                     'USERS_PASSWORDS_FLOWER': flower_users_passwords,
                     'JOB_MASTER': job_master,
                     'INSTANCE_TYPE': env.tags[env.host_string]['type']
                     }, use_sudo=True, use_jinja=True)

    # Copy env variable from postactivate
    run('python {source_dir}/deploy/cp_p2s.py -i {env_dir}/bin/postactivate -o {supervisor}'.format(
        source_dir=env.source_dir,
        env_dir=env.env_dir,
        supervisor=supervisor_file))

    # Simlink to conf file
    run_as_root('rm /etc/supervisor/supervisord.conf')
    run_as_root('ln -s {stack_dir}/{conf_dir}/supervisord.conf /etc/supervisor/supervisord.conf'.format(
        stack_dir=env.stack_dir,
        conf_dir=env.conf_dir,
    ))


@task
def restart_supervisor():
    # restart supervisor
    with settings(warn_only=True):
        pid = run('pgrep supervisor')
        if pid:
            run_as_root('kill -HUP {pid}'.format(pid=pid))
        else:   # just start it
            run('supervisord')


@task
def reload_nginx():
    run_as_root('sudo service nginx reload')


@task
def reload_supervisor():
    run_as_root('supervisorctl reread')
    run_as_root('supervisorctl update')


@task
def update_hosts_file(stack=STACK):
    conn = _create_connection(REGION)
    reservations = conn.get_all_instances(filters={'tag:stack': stack})

    context = {}
    for reservation in reservations:
        for instance in reservation.instances:
            if instance.public_dns_name:
                context[instance.public_dns_name] = {
                    'ip': instance.ip_address,
                    'private_ip': instance.private_ip_address,
                    'name': instance.tags.get('Name', ''),
                }
    # overwrite host_string ip to localhost
    context[env.host_string] = {
        'ip': '127.0.0.1',
        'private_ip': '127.0.0.1',
        'name': context[env.host_string]['name'],
    }

    # update host file if necessary
    reb = False
    for cont in context.values():
        if not cont['private_ip'] == '127.0.0.1':
            line_str = '{private_ip} {name} {name}.localdomain'.format(
                private_ip=cont['private_ip'],
                name=cont['name'],
            )
        else:
            line_str = '{private_ip} {name} {name}.localdomain localhost localhost.localdomain'.format(
                private_ip=cont['private_ip'],
                name=cont['name'],
            )

        if not files.contains('/etc/hosts', line_str):
            if files.contains('/etc/hosts', cont['private_ip']):  # update line if ip already in
                run_as_root("sed -i 's/.*{private_ip}.*/{new_line}/' /etc/hosts".format(
                    private_ip=cont['private_ip'],
                    new_line=line_str))
            elif files.contains('/etc/hosts', cont['name']):  # update line if name already in
                run_as_root("sed -i 's/.*{name}.*/{new_line}/' /etc/hosts".format(
                    name=cont['name'],
                    new_line=line_str))
            else:  # add line to top row
                run_as_root("sed -i -e '1i{new_line}\' /etc/hosts".format(
                    new_line=line_str))

    # update host_name
    if not files.contains('/etc/hostname', context[env.host_string]['name']):
        run_as_root('echo "{name}.localdomain" > /etc/hostname'.format(name=context[env.host_string]['name']))
        run_as_root('hostname {name}'.format(name=context[env.host_string]['name']))
        reb = True
    return reb

@task
def reboot_instance():
    reboot()

@task
def start_default():
    run('supervisorctl start celery-default')

@task
def start_consumers():
    run('supervisorctl start celery-consumers')

@task
def start_nlp():
    run('supervisorctl start celery-nlp')

@task
def start_ms():
    run('supervisorctl start celery-mostsimilar')

@task
def stop_default():
    run('supervisorctl stop celery-default')

@task
def stop_consumers():
    run('supervisorctl stop celery-consumers')

@task
def stop_nlp():
    run('supervisorctl stop celery-nlp')

@task
def stop_ms():
    run('supervisorctl stop celery-mostsimilar')

@task
def start_flower():
    run('supervisorctl start flower')

@task
def stop_flower():
    run('supervisorctl stop flower')

@task
def restart_flower():
    run('supervisorctl restart flower')

@task
def stop_all():
    run('supervisorctl stop all')

@task
def start_all():
    run('supervisorctl start all')


@task
def clear_nlp_data():
    run('rm -f /home/ubuntu/staging/source/nlp_data/mods/*')
    run('rm -f /home/ubuntu/staging/source/nlp_data/ms/*')
