#!/usr/bin/env python
"""
This script is executes when spot instances are going to be terminated by AWS
(e.g. due to price higher than bid). The script monitors
 http://169.254.169.254/latest/meta-data/spot/termination-time

If the termination signal is found, script triggers:
- cancel spot request
- start spot request with up to date AMI
"""

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import time
import requests
from distutils.version import LooseVersion
import os
from boto import utils
from boto.ec2 import connect_to_region, instance, connection
import re

DRY_RUN = True
SLEEP_TIME = 5  # sleep time between 2 checks
URL_CHECK = 'http://169.254.169.254/latest/meta-data/spot/termination-time'

if __name__ == '__main__':
    while True:
        time.sleep(5)
        resp = requests.get(URL_CHECK)
        if not resp.status_code == 404:  # instance is going to be terminated in 2min

            # Get instance data
            instance_id = utils.get_instance_identity()['document'][
                'instanceId']
            region = utils.get_instance_identity()['document']['region']

            # Get connection to AWS
            conn = connect_to_region(
                region_name=region,
                aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'),
                aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY')
            )

            # Retrieve Instance
            reservations = conn.get_all_instances(instance_ids=[instance_id])
            inst = reservations[0].instances[0]

            # Retrieve Spot request
            spot_id = inst.spot_instance_request_id
            spot = conn.get_all_spot_instance_requests(request_ids=[spot_id])[0]

            # Fetch some configuration
            security_groups = [grp.name for grp in inst.groups]

            # Get instance tags
            tags = inst.tags

            # Cancel spot request
            conn.cancel_spot_instance_requests(spot_id, dry_run=DRY_RUN)

            # Get corresponding AMI based on tags (+ latest version)
            filters = {
                'tag:stack': tags['stack'],
                'tag:layer': tags['layer'],
                'tag:role': tags['role'],
            }
            amis = conn.get_all_images(owners=['self'], filters=filters)
            # get last AMI version
            if len(amis) > 1:
                # get latest version
                versions = [ami.tags['version'] for ami in amis]
                last_v = LooseVersion(versions[0])
                for v in versions:
                    if LooseVersion(v) > last_v:
                        last_v = LooseVersion(v)
                last_version = last_v.vstring
                ami = \
                [ami for ami in amis if ami.tags['version'] == last_version][0]
            else:
                ami = amis[0]

            # Request spot instance
            new_spot = conn.request_spot_instances(spot.price,
                                                   ami.id,
                                                   count=1,
                                                   type=spot.type,
                                                   placement=inst.placement,
                                                   security_groups=security_groups,
                                                   instance_type=inst.instance_type,
                                                   key_name=inst.key_name,
                                                   ebs_optimized=True,
                                                   dry_run=DRY_RUN)[0]

            # Tag spot request
            new_spot.add_tags(ami.tags)
