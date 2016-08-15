#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import re
import os
import argparse
from subprocess import call

ROOT_DIR = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="bump version, merge to master, deploy")

    parser.add_argument("-v", "--version",
                        metavar='version',
                        help="Version number",
                        required=True,
                        type=str)
    parser.add_argument("-o", "--options",
                        help="Args to pass to deploy script",
                        default="pa",
                        type=str)
    args = parser.parse_args()

    # Create new release branch
    call(['git', 'checkout', '-b', 'release-{0}'.format(args.version)])

    # Bump version
    filename = os.path.join(ROOT_DIR, '__init__.py')
    with open(filename, "r") as f:
        lines = (line.rstrip() for line in f)
        altered_lines = []
        for line in lines:
            if re.match(r'__version__', line):
                altered_lines.append("__version__ = '{0}'".format(args.version))
            else:
                altered_lines.append(line)
    with open(filename, "w") as f:
        f.write('\n'.join(altered_lines) + '\n')
    call(['git', 'commit', '-am', 'bump to version {0}'.format(args.version)])

    # Merge with develop
    call(['git', 'checkout', 'develop'])
    call(['git', 'merge', '--no-ff', '--no-edit', 'release-{0}'.format(args.version)])

    # Merge with master
    call(['git', 'checkout', 'master'])
    call(['git', 'merge', '--no-ff', '--no-edit', 'release-{0}'.format(args.version)])

    # Push master
    call(['git', 'push'])

    # Deploy
    deploy_dir = './{0}'.format(os.path.join(ROOT_DIR, 'deploy'))
    call(['cd', deploy_dir])
    call(['workon', 'fab'])
    call(['./deploy.py', '-{0}'.format(args.options), 'production'])

    # Back to develop
    call(['cd', ROOT_DIR])
    call(['workon', 'paper'])
    call(['git', 'checkout', 'develop'])

