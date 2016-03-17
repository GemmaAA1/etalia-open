#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from boto import utils
from boto.ec2 import instance, connect_to_region

AWS_ACCESS_KEY_ID = 'AKIAJELT34R362VHB56A'
AWS_SECRET_ACCESS_KEY = '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'


def get_tags_from_spot_request():

    instance_id = utils.get_instance_identity()['document']['instanceId']
    region = utils.get_instance_identity()['document']['region']
    conn = connect_to_region(
        region_name=region,
        aws_access_key_id=AWS_ACCESS_KEY_ID,
        aws_secret_access_key=AWS_SECRET_ACCESS_KEY)
    inst = instance.Instance(connection=conn)
    inst.id = instance_id
    inst.update()
    spot_id = inst.spot_instance_request_id
    tags = conn.get_all_tags(filters={'resource-type': 'spot-instances-request', 'resource-id': spot_id})
    # remove Name tag
    tags.pop('Name')
    for tag in tags:
        inst.add_tag(tag.name, tag.value)

    # generate unique name


if __name__ == '__main__':
    get_tags_from_spot_request()
