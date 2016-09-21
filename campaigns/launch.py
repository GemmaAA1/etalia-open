#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import sys
import argparse
import django
import json


def fetch_path():
    # search for config/ is parent tree directory and setup django
    cpath = os.path.abspath(__file__)
    while cpath and not os.path.isdir('config'):
        os.chdir('..')
    sys.path.append(os.path.abspath(os.getcwd()))
    return os.path.abspath(os.getcwd())


def check_config(config):

    # check for duplicate names
    names = [c.get('name') for c in config]
    assert len(set(names)) == len(names)


if __name__ == '__main__':

    fetch_path()
    django.setup()
    from campaigns.utils import get_or_create_campaign
    from etalia.core.emails import Email
    from django.conf import settings

    # Input
    parser = argparse.ArgumentParser(description='Launch a mailing campaign\n')
    parser.add_argument("-c",
                        help="path to the inline email template file (default: campaigns.json)",
                        metavar='file_path',
                        dest='config',
                        type=str,
                        required=False,
                        default=str(settings.ROOT_DIR('campaigns/campaigns.json')))
    parser.add_argument("-n",
                        help="Name of the campaign as defined in the config file",
                        metavar='name',
                        dest='name',
                        type=str,
                        required=True)

    args = parser.parse_args()

    # Load config file
    with open(args.config, 'r') as config_file:
        config = json.load(config_file)

    # Check integrity of campaign config
    check_config(config)

    # Check if name in config
    assert args.name in [c.get('name') for c in config]

    # Get campaign
    campaign = [c for c in config if c.get('name') == args.name][0]

    # Create campaign
    mg_campaign = get_or_create_campaign(campaign.get('name'))

    # Parse emails list
    with open(str(settings.ROOT_DIR(campaign.get('email_list'))), 'r') as file:
        emails = file.readlines()
        emails = [email.rstrip() for email in emails]

    # Instantiate email class
    email = Email(
        template=str(settings.ROOT_DIR(campaign.get('email_template'))),
        cids=campaign.get('cids', {}),
        tags=campaign.get('tags', []),
        metadata=campaign.get('metadata', {}),
        subject=campaign.get('subject', ''),
        from_email=campaign.get('from_email'),
        to=[],
        reply_to=campaign.get('reply_to', []),
        extra_ctx=campaign.get('extra_ctx', {}))

    # send email
    for email in emails:
        email.to = [email]
        email.send()


