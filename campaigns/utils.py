# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import requests

from django.conf import settings

MAILGUN_API_KEY = settings.ANYMAIL['CUSTOMMAILGUN_API_KEY']
MAILGUN_API_URL = 'https://api.mailgun.net/v3/{domain}/campaigns'.format(
    domain=settings.ANYMAIL['CUSTOMMAILGUN_SENDER_DOMAIN'])


def get_or_create_campaign(name, campaign_id=None):

    # Check campaign based on id
    if campaign_id:
        resp = get_campaign(campaign_id)
        if resp.status_code == 200:
            if resp.json().name == name:
                return resp.json()
            else:
                raise ValueError('name and campaign_id does not match')

    # Check campaign based on name
    data = get_campaigns().json()
    for c in data['items']:
        if name == c.get('name'):
            if campaign_id == c.get('id'):
                return c
            elif campaign_id is None:
                return c

    # If not found create
    resp = create_campaign(name, campaign_id=campaign_id)
    if resp.status_code == 200:
        return resp.json().get('campaign')

    return None


def get_campaigns():
    """Return list of campaigns"""
    return requests.get(
        MAILGUN_API_URL,
        auth=('api', MAILGUN_API_KEY)
    )


def get_campaign(campaign_id):
    """Return campaign"""
    return requests.get(
        '{base}/{id}'.format(base=MAILGUN_API_URL, id=campaign_id),
        auth=('api', MAILGUN_API_KEY)
    )


def create_campaign(name, campaign_id=None):
    """Create new campaign"""
    data = {'name': name}
    if campaign_id:
        data['id'] = campaign_id
    return requests.post(MAILGUN_API_URL, data, auth=('api', MAILGUN_API_KEY))


def delete_campaign(campaign_id):
    """Delete campaign"""
    return requests.delete('{base}/{id}'.format(base=MAILGUN_API_URL,
                                                id=campaign_id))


def get_events_history(campaign_id, **kwargs):
    """Retrieve events history for a campaign

    Kwargs:
        event (str): Name of event to filter by - clicked, opened,
            unsubscribed, bounced, complained. (optional)
        recipient (str): Address of recipient to filter by. (optional)
        country (str): Two-letters country code to filter by. For example: US,
            RU. (optional)
        region (str): Name or code of region to filter by. For example, US
            states: CA, OH. (optional)
        limit (int): Maximum number of records to return (100 at most).
            (optional)
        page (int): Number of page to retrieve. (optional)
        count (int): Toggles whether to return actual data or just count of
            records. (optional)
    """
    return requests.get(
        '{base}/{id}/events'.format(
            base=MAILGUN_API_URL,
            id=campaign_id,),
        auth=('api', MAILGUN_API_KEY),
        params=kwargs
    )


def get_campaign_stats(campaign_id, **kwargs):
    """Retrieve general stats for a campaign

    Kwargs:
        groupby (str): Grouping parameter - domain or daily_hour. (optional)
    """
    return requests.get(
        '{base}/{id}/stats'.format(base=MAILGUN_API_URL, id=campaign_id),
        auth=('api', MAILGUN_API_KEY),
        params=kwargs
    )
