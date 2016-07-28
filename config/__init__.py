import boto3
from .celery import celery_app


def get_dns_name_based_on_role(role):

    ec2 = boto3.resource('ec2')
    instances = list(ec2.instances.all())
    props = [(i.private_dns_name, i.tags) for i in instances]
    for prop in props:
        for tag in prop[1]:
            if tag['Key'] == 'role' and role in tag['Value']:
                return prop[0]
