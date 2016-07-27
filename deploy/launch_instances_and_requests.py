#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
import os
import json
from aws import connect_ec2, get_latest_ami, get_tag_val, get_instance_name

env = environ.Env()
CONFIG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'config.json')
DRY_RUN = False
SPOT_MAX_ACTIVE_TIME_WAIT = 30  # sec


if __name__ == '__main__':

    # Load configuration
    with open(CONFIG_FILE) as file:
        conf = json.load(file)

    # Connect to EC2
    ec2 = connect_ec2()

    for r in conf:
        # Create instance
        if r.get('Instance'):
            props = r.get('Instance')
            props['DryRun'] = True if DRY_RUN else False
            tags = r.get('Tags')
            ami = get_latest_ami(ec2, r.get('Tags'))[0]
            ins = ec2.create_instances(ImageId=ami.id, **props)[0]

            # Tag instance
            tags.append({"Key": "version", "Value": get_tag_val(ami, 'version')})
            tags.append({"Key": "Name", "Value": get_instance_name(ec2, tags)})
            # tags = ins.create_tags(Tags=tags)

        elif r.get('Spot'):
            props = r.get('Spot')
            props['DryRun'] = True if DRY_RUN else False
            tags = r.get('Tags')
            ami = get_latest_ami(ec2, tags)[0]
            props['LaunchSpecification']['ImageId'] = ami.id
            spot = ec2.meta.client.request_spot_instances(**props)
