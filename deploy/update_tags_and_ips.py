#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Update instance tags either based on AMI tag (if tags is empty, i.e. this is
a newly launched instance) or based on instance current tags and newly fetched
data (i.e. after a new deploy for instance)"""

from __future__ import unicode_literals, absolute_import
from aws import connect_ec2, get_local_instance_id, tag_instance, \
    associate_elastic_ip
from time import sleep

if __name__ == '__main__':

    # Get local instance id
    instance_id = get_local_instance_id()

    # Connect to EC2
    ec2 = connect_ec2()

    # Build tags
    tag_instance(ec2, instance_id)

    # Associate elastic IP
    associate_elastic_ip(ec2, instance_id)

    # Go to sleep for 10 sec
    sleep(20)




