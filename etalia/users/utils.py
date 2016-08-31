# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import datetime
from django.conf import settings
from django.template.loader import get_template
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from django.contrib.auth import get_user_model
from django.utils import timezone

from .models import UserInvited

User = get_user_model()


def send_invite_email(email_to=None,
                      root_url=None,
                      on_behalf=None):
    subject = 'An invitation to try Etalia'
    to = [email_to]
    from_email =    'nicolas.pannetier@etalia.io'
    if on_behalf:
        user_name= '{first} {last}'.format(
                first=on_behalf.first_name,
                last=on_behalf.last_name
        )
    else:
        user_name = None
    ctx = {
        'bucket_url': settings.EMAIL_STATIC_BUCKET,
        'root_url': root_url,
        'on_behalf': user_name,
    }
    text_content = ''
    html_content = get_template(settings.INVITE_EMAIL_TEMPLATE)\
        .render(Context(ctx))
    email = EmailMultiAlternatives(subject, text_content, to=to, from_email=from_email)
    email.attach_alternative(html_content, "text/html")
    # email.content_subtype = 'html'
    email.send()

    # save to database
    UserInvited.objects.create(from_user=on_behalf, to_email=email_to)


def send_periodic_recommendation_email(user_id):

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
        subject = 'Recommendations from Etalia'
        from_email = 'etalia@etalia.io'
        to = [user.email]

        ctx = {
            'bucket_url': settings.EMAIL_STATIC_BUCKET,
            'papers': papers,
            'root_url': 'http://etalia.io',
        }

        text_content = ''
        html_content = get_template(settings.PERIODIC_RECOMMENDATION_TEMPLATE)\
            .render(Context(ctx))
        email = EmailMultiAlternatives(subject,
                                       text_content,
                                       to=to,
                                       from_email=from_email)
        email.attach_alternative(html_content, "text/html")
        email.send()


def next_weekday(date, weekday):
    days_ahead = weekday - date.weekday()
    if days_ahead <= 0:   # Target day already happened this week
        days_ahead += 7
    return date + datetime.timedelta(days_ahead)
