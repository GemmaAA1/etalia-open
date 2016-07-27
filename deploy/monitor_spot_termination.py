#!/usr/bin/env python
"""
This script is executes when spot instances are going to be terminated by AWS
(e.g. due to price higher than bid). The script monitors
 http://169.254.169.254/latest/meta-data/spot/termination-time

If the termination signal is found, script triggers:
- cancel spot request
- start a new spot request with latest AMI
"""

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import time
import requests
from aws import connect_ec2, get_local_instance_id, get_latest_ami, \
    URL_SPOT_TERMINATION_CHECK

DRY_RUN = False
SLEEP_TIME = 5  # in s


def going_to_termination():
    """Check for termination signal from aws spot instance"""
    resp = requests.get(URL_SPOT_TERMINATION_CHECK)
    if resp.status_code == 404:
        return False
    else:
        return True


if __name__ == '__main__':

    while True:

        time.sleep(5)

        if going_to_termination():

            # connect
            ec2 = connect_ec2()

            # Get local instance
            instance_id = get_local_instance_id()
            instance = list(ec2.instances.filter(InstanceIds=[instance_id]))[0]
            spot_id = instance.spot_instance_request_id

            # Describe spot instance
            spot = ec2.meta.client.describe_spot_instance_requests(
                SpotInstanceRequestIds=[spot_id])

            # Cancel related spot instance
            ec2.meta.client.cancel_spot_instance_requests(
                SpotInstanceRequestIds=[spot_id])

            # start a spot request with updated ami
            ami = get_latest_ami(ec2, instance.tags)[0]
            props = spot['SpotInstanceRequests'][0]
            props['LaunchSpecification']['ImageId'] = ami.id
            ec2.meta.client.request_spot_instances(**props)

            # terminate instance
            ec2.meta.client.terminate_instances(InstanceIds=[instance_id, ])
