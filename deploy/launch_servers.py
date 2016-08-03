#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import environ
import os
import json
from aws import connect_ec2, get_latest_ami, get_tag_val, get_instance_name, \
    dict2tags, tags2dict

env = environ.Env()
# CONFIG_FILE = 'config.json'
CONFIG_FILE = 'config_fleet.json'
DRY_RUN = True
SPOT_MAX_ACTIVE_TIME_WAIT = 30  # sec


def launch(k, config=None):

    # Load configuration
    if not config:
        with open(CONFIG_FILE) as file:
            config = json.load(file)

    instance_type = config[k]['type']
    props = config[k]['props']
    ami_tags = config[k]['ami-tags']

    # overwrite DryRung prop
    props['DryRun'] = True if DRY_RUN else False

    # Connect to EC2
    ec2 = connect_ec2()

    # Get ami
    ami = get_latest_ami(ec2, dict2tags(tags2dict(ami_tags), filters=True))[0]

    if instance_type == 'instance':
        launch_instance(props, ami.id)
    elif instance_type == 'spot':
        launch_spot(props, ami.id)
    elif instance_type == 'fleet':
        launch_fleet(props, ami.id)
    else:
        raise ValueError('unknown instance type {0}'.format(instance_type))


def launch_instance(props, ami_id):
    # Connect to EC2
    ec2 = connect_ec2()
    # Create instance
    ins = ec2.create_instances(ImageId=ami_id, **props)[0]
    return ins


def launch_spot(props, ami_id):
    # Connect to EC2
    ec2 = connect_ec2()
    ['LaunchSpecification']['ImageId'] = ami_id
    # Request spot instance
    spot = ec2.meta.client.request_spot_instances(**props)
    return spot


def launch_fleet(props, ami_id):
    # Connect to EC2
    ec2 = connect_ec2()

    instance_types = props['SpotFleetRequestConfig'].pop('InstanceTypes')
    key_name = props['SpotFleetRequestConfig'].pop('KeyName')
    security_groups = props['SpotFleetRequestConfig'].pop("SecurityGroups")

    props['SpotFleetRequestConfig']['LaunchSpecifications'] = []
    for itype in instance_types:
        props['SpotFleetRequestConfig']['LaunchSpecifications'].append(
            {
                "ImageId": ami_id,
                "InstanceType": itype,
                "KeyName": key_name,
                "SecurityGroups": security_groups
            }
        )

    # Request spot instance
    spot = ec2.meta.client.request_spot_fleet(**props)
    return spot


if __name__ == '__main__':

    # Load configuration
    with open(CONFIG_FILE) as file:
        conf = json.load(file)

    for key, val in conf.items():
        launch(key, config=conf)
