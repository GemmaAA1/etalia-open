# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import argparse
from etalia.users.utils import send_invite_email
from time import sleep


email_file = 'beta_emails.txt'


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Send invite email in batch')
    parser.add_argument("-f", "--file",
                        help="email file list",
                        action="store_true")
    parser.add_argument("-b", "--behalf",
                        help="email file list",
                        action="store_true",
                        default="etalia@etalia.io")

    args = parser.parse_args()

    # Get email list
    with open(args.file, 'r') as f:
        emails = f.readlines()
    emails = [email.strip() for email in emails]

    # Get behalf
    on_behalf = args.behalf

    for email in emails:
        sleep(0.5)
        send_invite_email(
            email_to=email,
            on_behalf=on_behalf,
            root_url='etalia.io'
        )
