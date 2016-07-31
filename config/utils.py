# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import boto3


def get_dns_name_based_on_role(role):

    ec2 = boto3.resource('ec2')
    instances = list(ec2.instances.all())
    props = [(i.private_dns_name, i.tags) for i in instances]
    for prop in props:
        for tag in prop[1]:
            if tag['Key'] == 'role' and role in tag['Value']:
                return prop[0]


def get_private_ip_based_on_role(role):

    ec2 = boto3.resource('ec2')
    instances = list(ec2.instances.all())
    props = [(i.private_ip_address, i.tags) for i in instances]
    for prop in props:
        for tag in prop[1]:
            if tag['Key'] == 'role' and role in tag['Value']:
                return prop[0]