# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import datetime
from django.conf import settings
from django.contrib.auth import get_user_model
from django.utils import timezone

from etalia.core.emails import Email
from .models import UserInvited

User = get_user_model()


def send_invite_email(email_to=None, on_behalf=None):
    """Send invitation email"""

    # Get data
    if on_behalf:
        user_name = '{first} {last}'.format(first=on_behalf.first_name,
                                            last=on_behalf.last_name)
    else:
        user_name = None

    # Instantiate
    email = Email(
        template=settings.INVITE_EMAIL_TEMPLATE,
        cids={'logo_cid': 'beta_etalia_logo.png',
              'book_cid': 'book.png',
              'engage_cid': 'engage.png',
              'target_cid': 'target.png'},
        tags=['invite', 'peer-to-peer'],
        metadata={'from_user': on_behalf.id},
        subject='An invitation to try Etalia',
        from_email='nicolas.pannetier@etalia.io',
        to=[email_to],
        reply_to=['contact@etalia.io'],
        extra_ctx={'on_behalf': user_name})
    # send email
    email.send()

    # save to database
    UserInvited.objects.create(from_user=on_behalf, to_email=email_to)


def send_periodic_recommendation_email(user_id):
    """Send etalia digest email"""

    # Get data
    users = User.objects.filter(id=user_id).select_related('settings')
    user = users[0]
    date_7d_before = timezone.now() - datetime.timedelta(days=7)
    if hasattr(user, 'userperiodicemail') and user.userperiodicemail.last_sent_on:
        date_since = user.userperiodicemail.last_sent_on
        if date_since > date_7d_before:
            date_since = date_7d_before
    else:
        date_since = datetime.datetime.now() - \
                     datetime.timedelta(days=user.settings.email_digest_frequency)
    papers = user.streams.first().papers\
        .filter(streampapers__date__gte=date_since)\
        .order_by('-streampapers__score')
    papers = papers.select_related('altmetric')
    papers = papers.select_related('journal')
    papers = papers.prefetch_related('authors')
    papers = list(papers[:settings.PERIODIC_RECOMMENDATION_NUMBER_PAPERS])

    if papers:

        email = Email(
            template=settings.PERIODIC_RECOMMENDATION_TEMPLATE,
            cids={'logo_cid': 'etalia_digest.png'},
            tags=['digest'],
            metadata={'user': user.id},
            subject='Recommendations from Etalia',
            from_email='etalia@etalia.io',
            to=[user.email],
            reply_to=['contact@etalia.io'],
            extra_ctx={'papers': papers})

        email.send()


def next_weekday(date, weekday):
    days_ahead = weekday - date.weekday()
    if days_ahead <= 0:   # Target day already happened this week
        days_ahead += 7
    return date + datetime.timedelta(days_ahead)
