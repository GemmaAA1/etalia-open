#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from boto import utils
from boto.ec2 import instance, connect_to_region

AWS_ACCESS_KEY_ID = 'AKIAJELT34R362VHB56A'
AWS_SECRET_ACCESS_KEY = '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'


def get_tags_from_spot_request():

    # Get instance id, region
    instance_id = utils.get_instance_identity()['document']['instanceId']
    region = utils.get_instance_identity()['document']['region']

    # Connect to AWS
    conn = connect_to_region(
        region_name=region,
        aws_access_key_id=AWS_ACCESS_KEY_ID,
        aws_secret_access_key=AWS_SECRET_ACCESS_KEY)

    # Get instance object
    inst = instance.Instance(connection=conn)
    inst.id = instance_id
    inst.update()

    # Get related spot instance
    spot_id = inst.spot_instance_request_id
    # get tags from spot
    tags = conn.get_all_tags(filters={'resource-type': 'spot-instances-request', 'resource-id': spot_id})

    # Add tags to instance
    tags.pop('Name')
    inst.add_tags(tags)

    # Generate unique name
    # grab base name from ami
    ami = conn.get_all_images(image_ids=[inst.image_id])[0]
    base_name = ami.name
    # fetch all current name
    rs = conn.get_all_reservations()
    inst_names = [r.instances[0].tags['Name'] for r in rs]
    # Compare names and increment
    count = 0
    inst_name = '{base}.id{id:03d}'.format(base=base_name, id=count)
    while inst_name in inst_names:
        count += 1
        inst_name = '{base}.id{id:03d}'.format(base=base_name, id=count)
    inst.add_tag({'Name': inst_name})


if __name__ == '__main__':
    get_tags_from_spot_request()
