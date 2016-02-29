#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import os
import environ
from subprocess import call
from fabric_slack_tools import *

SLACK_WEB_HOOK = "https://hooks.slack.com/services/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"


def send_deploy_version_message():
    # init slack web hook
    init_slack(SLACK_WEB_HOOK)
    ROOT_DIR = environ.Path(__file__) - 2  # (/a/myfile.py - 2 = /)
    # Get app version from root __init__
    version = get_version(str(ROOT_DIR.path()))
    send_slack_message("Deploying Etalia {0}...".format(version),
                       channel="#general",
                       username="deployment-bot",
                       web_hook_url=SLACK_WEB_HOOK)


def get_version(package):
    """
    Return package version as listed in `__version__` in `init.py`.
    """
    init_py = open(os.path.join(package, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)


def checkout_master_and_push_recompiled_assets():
    # Git checkout master
    call(["git", "checkout", "master"])
    # Recompiled assets and push to repo
    call('gulp')
    call(['git', 'commit', '-am', 'recompiled assets prior deployment'])
    call(['git', 'push'])


def deploy(stack='production'):

    send_deploy_version_message()
    checkout_master_and_push_recompiled_assets()
    call(['fab', 'set_hosts:{stack},*,*'.format(stack=stack), '-P',  'deploy'])

if __name__ == '__main__':
    deploy()
