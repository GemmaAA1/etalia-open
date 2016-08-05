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
import subprocess
from aws import connect_ec2, get_local_instance_id, get_latest_ami, \
    URL_SPOT_TERMINATION_CHECK, tags2dict, dict2tags
import datetime

DRY_RUN = False
SLEEP_TIME = 5  # in s
VALID_UNTIL = datetime.datetime(2020, 1, 1)


SPOT_TEMPLATE = {
    'DryRun': False,
    'SpotPrice': '0.1',
    'InstanceCount': 1,
    'Type': 'one-time',
    'AvailabilityZoneGroup': 'us-west-2b',
    'LaunchSpecification': {
        'KeyName': 'npannetier-key-pair',
        'SecurityGroups': ['paperstream_sg'],
        'InstanceType': 'm4.xlarge',
        'Placement': {
            'AvailabilityZone': 'us-west-2b'
        }
    }
}


def going_to_termination():
    """Check for termination signal from aws spot instance"""
    resp = requests.get(URL_SPOT_TERMINATION_CHECK)
    if resp.status_code == 404:
        return False
    else:
        return True


def get_fleet_id_from_instance_id(instance_id):

    # connect
    ec2 = connect_ec2()

    # Get fleet request
    sfrcs = ec2.meta.client.describe_spot_fleet_requests()["SpotFleetRequestConfigs"]

    # Filter active fleet request:
    active_ids = []
    for sfrc in sfrcs:
        if sfrc['SpotFleetRequestState'] == 'active':
            active_ids.append(sfrc['SpotFleetRequestId'])

    # Loop through active fleet to spot instance_id
    for sfid in active_ids:
        res = ec2.meta.client.describe_spot_fleet_instances(
            SpotFleetRequestId=sfid
        )['ActiveInstances']
        for r in res:
            if instance_id in r['InstanceId']:
                return sfid

    raise None


def cancel_and_restart():

    # connect
    ec2 = connect_ec2()

    # Get local instance
    instance_id = get_local_instance_id()
    instance = list(ec2.instances.filter(InstanceIds=[instance_id]))[0]

    # Get relate spot and fleet ids
    spot_id = instance.spot_instance_request_id
    fleet_id = get_fleet_id_from_instance_id(instance_id)

    # Get latest AMI for instance
    d = tags2dict(instance.tags)
    d.pop('Name', None)
    for key in list(d.keys()):
        if key.startswith('aws:'):
            d.pop(key)
    ami = get_latest_ami(ec2, dict2tags(d, filters=True))[0]

    if fleet_id:  # spot instance from fleet

        # Get fleet config
        fleet = ec2.meta.client.describe_spot_fleet_requests(
            SpotFleetRequestIds=[fleet_id]
        )

        # Cancel related fleet request
        ec2.meta.client.cancel_spot_fleet_requests(
            SpotFleetRequestIds=[fleet_id],
            TerminateInstances=False
        )

        # Update config with new amis
        props = fleet['SpotFleetRequestConfigs'][0]['SpotFleetRequestConfig']
        for ins in props['LaunchSpecifications']:
            ins['ImageId'] = ami.id
            ins.pop('BlockDeviceMappings', None)

        props.pop('ValidFrom', None)
        props.pop('ClientToken', None)
        props['ValidUntil'] = VALID_UNTIL

        # Start a new fleet with latest ami
        ec2.meta.client.request_spot_fleet(DryRun=DRY_RUN,
                                           SpotFleetRequestConfig=props)

    elif spot_id:  # spot instance not from fleet
        # Describe spot instance
        spot = ec2.meta.client.describe_spot_instance_requests(
            SpotInstanceRequestIds=[spot_id])

        # Cancel spot request instance
        ec2.meta.client.cancel_spot_instance_requests(
            SpotInstanceRequestIds=[spot_id])

        # start a spot request with latest ami
        props = SPOT_TEMPLATE
        last_spot = spot['SpotInstanceRequests'][0]
        props['LaunchSpecification']['ImageId'] = ami.id
        props['SpotPrice'] = last_spot['SpotPrice']
        props['LaunchSpecification']['InstanceType'] = \
            last_spot['LaunchSpecification']['InstanceType']
        props['DryRun'] = DRY_RUN
        ec2.meta.client.request_spot_instances(**props)
    else:
        raise EnvironmentError('Not a spot instance')

    # Try a clean service shutdown
    subprocess.call(["supervisorctl", "stop", "all"])


if __name__ == '__main__':

    while True:

        time.sleep(5)

        if going_to_termination():

            cancel_and_restart()

            break
