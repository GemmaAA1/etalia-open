# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import boto3
import re
import os
import environ
import requests
from distutils.version import LooseVersion

INSTANCE_TAGS_KEY = ['Name', 'stack', 'layer', 'role', 'version']
INSTANCE_NAME_PATTERN = ['stack', 'layer', 'role', 'version']
IMAGE_NAME_PATTERN = ['version', 'stack', 'layer', 'role']
URL_INSTANCE_ID_CHECK = 'http://169.254.169.254/latest/meta-data/instance-id'
URL_SPOT_TERMINATION_CHECK = \
    'http://169.254.169.254/latest/meta-data/spot/termination-time'
ROOT_DIR = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def connect_ec2():
    """Return ec2 resource"""
    return boto3.resource('ec2')


def get_local_instance_id():
    response = requests.get(URL_INSTANCE_ID_CHECK)
    return response.text


def get_latest_ami(ec2, filters):
    """Return latest version of AMIs matching tags"""

    # Retrieve AMIs
    amis = list(ec2.images.filter(Filters=filters, Owners=['self']))

    # Filter AMIs based on latest version tag
    try:
        if len(amis) > 1:
            # extract version from tags
            versions = []
            for ami in amis:
                versions.append(get_tag_val(ami, 'version'))

            # get latest version
            last_v = LooseVersion(versions[0])
            for v in versions:
                if LooseVersion(v) > last_v:
                    last_v = LooseVersion(v)
            last_version = last_v.vstring

            # filter latest version
            latest_amis = []
            for i, v in enumerate(versions):
                if v == last_version:
                    latest_amis.append(amis[i])
            return latest_amis
        else:
            return [amis[0]]
    except IndexError:
        raise 'No AMI found with tags {0}'.format(filters)


def get_tag_val(obj, key):
    """Return ami version found in tags"""
    for tag in obj.tags:
        if tag.get('Key') == key:
            return tag.get('Value')
    return 'unknown'


def get_all_running_instances(ec2):
    """Return list of running instances"""
    instances = ec2.instances.filter(Filters=[{'Name': 'instance-state-name',
                                               'Values': ['running']}])
    return list(instances)


def get_instance_name(ec2, tags):
    """Build instance name as stack-layer-role-version-#"""

    # Build base name
    name_list = []
    for a in INSTANCE_NAME_PATTERN:
        for tag in tags:
            if tag.get('Key') == a:
                name_list.append(tag.get('Value'))
    base = '-'.join(name_list)
    base = base.replace('.', '-')
    base = base.replace('/', '-')

    # fetch all current name
    instances = get_all_running_instances(ec2)
    inst_names = [get_tag_val(i, 'Name') for i in instances if not 'unknow']

    # Compare names and increment
    count = 0
    inst_name = '{base}--{id:03d}'.format(base=base, id=count)
    while inst_name in inst_names:
        count += 1
        inst_name = '{base}--{id:03d}'.format(base=base, id=count)

    return inst_name


def get_image_name(tags):
    """Build instance name as stack-layer-role-version"""

    # Build base name
    name_list = []
    for a in IMAGE_NAME_PATTERN:
        for tag in tags:
            if tag.get('Key') == a:
                name_list.append(tag.get('Value'))
    base = '_'.join(name_list)
    base = base.replace('.', '-')
    base = base.replace('/', '-')

    return base


def clean_tags(tags, drop=None):
    """Return list of tags compatible with keys in INSTANCE_TAGS_KEY but drop"""
    if drop is not None and not isinstance(drop, list):
        raise ValueError('drop kwargs must be a list of str')
    if drop is None:
        drop = []
    valid_keys = [key for key in INSTANCE_TAGS_KEY if key not in drop]
    return [tag for tag in tags if tag['Key'] in valid_keys]


def tag_instance(ec2, instance_id):
    """Tag instance"""
    instance = list(ec2.instances.filter(InstanceIds=[instance_id]))[0]
    try:
        ami = list(ec2.images.filter(ImageIds=[instance.image_id]))[0]
        tags = ami.tags
    except IndexError:  # AMI id has been removed, use local tags
        tags = instance.tags
    tags = clean_tags(tags, drop=['Name', 'version'])
    tags.append({"Key": "version", "Value": get_etalia_version()})
    tags.append({"Key": "Name", "Value": get_instance_name(ec2, tags)})

    # Update tags
    tags = instance.create_tags(Tags=tags)

    return tags


def dict2tags(dic, filters=False):
    if filters:
        return list(map(lambda x: {'Name': 'tag:' + x[0], 'Values': [x[1]]}, list(dic.items())))
    else:
        return map(lambda x: {'Key': x[0], 'Value': x[1]}, list(dic.items()))


def tags2dict(tags):
    return dict(map(lambda x: (x['Key'], x['Value']), tags or []))


def get_etalia_version():
    init_py = open(os.path.join(ROOT_DIR, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)


def associate_elastic_ip(ec2, instance_id):
    """Associate elastic based on instance tag"""

    from .elastic_ip import ELASTIC_IP_MAPPING

    instance = list(ec2.instances.filter(InstanceIds=[instance_id]))[0]
    tags = tags2dict(instance.tags)
    roles = tags.get('role', [])

    for k, props in ELASTIC_IP_MAPPING.items():
        if k in roles:
            rsp = ec2.meta.client.associate_address(
                DryRun=False,
                InstanceId=instance_id,
                AllocationId=props.get('allocation_id'),
                AllowReassociation=True
            )


def get_fleet_request_id_from_instance_id(instance_id):

    ec2 = connect_ec2()

    # Get all active fleet requests
    resp = ec2.meta.client.describe_spot_fleet_requests()
    fleet_configs = resp['SpotFleetRequestConfigs']

