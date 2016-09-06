#!/usr/bin/env python
from __future__ import unicode_literals, absolute_import
import sys
import os
import argparse
from time import sleep


def fetch_path():
    # search for config/ is parent tree directory and setup django
    cpath = os.path.abspath(__file__)
    while cpath and not os.path.isdir('config'):
        os.chdir('..')
    sys.path.append(os.path.abspath(os.getcwd()))

    return os.path.abspath(os.getcwd())


if __name__ == '__main__':

    # add project root to python path
    root_path = fetch_path()

    import django
    django.setup()
    from etalia.users.utils import send_invite_email
    from django.contrib.auth import get_user_model
    User = get_user_model()

    parser = argparse.ArgumentParser(description='Send invite email in batch')
    parser.add_argument("-f", "--file",
                        help="email file list",
                        type=str)
    parser.add_argument("-b", "--behalf",
                        help="register user email",
                        type=str,
                        default="etalia@etalia.io")

    args = parser.parse_args()

    # Get email list
    file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'scripts/routines',
                        args.file)
    print('loading file: {}'.format(file))
    with open(file, 'r') as f:
        emails = f.readlines()
    emails = [email.strip() for email in emails]

    # Get behalf
    user = User.objects.get(email=args.behalf)

    for email in emails:
        sleep(0.5)
        send_invite_email(
            email_to=email,
            on_behalf=user,
            root_url='https://etalia.io'
        )
