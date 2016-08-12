# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from etalia.users.utils import send_invite_email
from time import sleep


email_file = 'beta_emails.txt'


if __name__ == '__main__':
    with open(email_file, 'r') as f:
        emails = f.readlines()

    emails = [email.strip() for email in emails]

    for email in emails:
        sleep(0.5)
        send_invite_email(
            email_to=email,
            on_behalf=None,
            root_url='etalia.io'
        )
