#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
To deploy etalia on AWS.

>> ./deploy.py

Requirements:
- python 2.7    (fabric is python 2.7 only)
- pip install -r deploy_requirements.txt

"""

from __future__ import unicode_literals, absolute_import

import re
import os
import environ
import argparse
from subprocess import call
from time import sleep
from fabric_slack_tools import *

SLACK_WEB_HOOK = "https://hooks.slack.com/services/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"
init_slack(SLACK_WEB_HOOK)


def send_deploy_version_message(stack, done=False):
    ROOT_DIR = environ.Path(__file__) - 2  # (/a/myfile.py - 2 = /)
    # Get app version from root __init__
    version = get_version(str(ROOT_DIR.path()))
    if done:
        mess = "Deploying Etalia {vers} on {stack} stack is done".format(vers=version, stack=stack)
    else:
        mess = "Deploying Etalia {vers} on {stack} stack...".format(vers=version, stack=stack)

    send_slack_message(mess,
                       channel="#general",
                       username="deployment-bot",
                       web_hook_url=SLACK_WEB_HOOK)

def get_version(package):
    """
    Return package version as listed in `__version__` in `init.py`.
    """
    init_py = open(os.path.join(package, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)


def checkout_master():
    # Git checkout master
    call(["git", "checkout", "master"])


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("stack", help="(str) stack name to be deploy", type=str)
    parser.add_argument("-p", "--parallel", help="deploy in parallel", action="store_true")
    parser.add_argument("-s", "--slack", help="verbose to Slack", action="store_true")
    parser.add_argument("-a", "--amis", help="create AMIs", action="store_true")
    args = parser.parse_args()

    if args.slack:
        send_deploy_version_message(args.stack)

    if args.parallel:
        # Go on maintenance
        call(['fab', 'set_hosts:{stack},apps,*'.format(stack=args.stack), '-P', 'go_on_maintenance'])
        # parallel deploy
        if args.slack:
            call(['fab', 'set_hosts:{stack},*,*'.format(stack=args.stack), '-P', 'deploy_verbose'])
        else:
            call(['fab', 'set_hosts:{stack},*,*'.format(stack=args.stack), '-P', 'deploy'])
    else:
        # Go on maintenance
        call(['fab', 'set_hosts:{stack},apps,*'.format(stack=args.stack), 'go_on_maintenance'])
        # parallel deploy
        if args.slack:
            call(['fab', 'set_hosts:{stack},*,*'.format(stack=args.stack), 'deploy_verbose'])
        else:
            call(['fab', 'set_hosts:{stack},*,*'.format(stack=args.stack), 'deploy'])

    # Create AMIs
    if args.amis:
        print('Let'' wait few seconds before AMI registration...')
        sleep(20)    # for instance tag to be updated
        ROOT_DIR = environ.Path(__file__) - 2  # (/a/myfile.py - 2 = /)
        version = get_version(str(ROOT_DIR.path()))
        if args.slack:
            send_slack_message('Registering AMIs {v}'.format(v=version),
                       channel="#general",
                       username="deployment-bot",
                       web_hook_url=SLACK_WEB_HOOK)
        call(['fab', 'create_amis:{stack}'.format(stack=args.stack)])

    # Warm-up dispatchers
    print('Let'' wait few seconds before warming up dispatcher...')
    sleep(15)
    call(['fab', 'set_hosts:{stack},jobs,*'.format(stack=args.stack), 'warm_up_dispatchers'])

    if args.parallel:
        # Go off maintenance
        call(['fab', 'set_hosts:{stack},apps,*'.format(stack=args.stack), '-P', 'go_off_maintenance'])
    else:
        # Go off maintenance
        call(['fab', 'set_hosts:{stack},apps,*'.format(stack=args.stack), 'go_off_maintenance'])

    if args.slack:
        send_deploy_version_message(args.stack, done=True)
