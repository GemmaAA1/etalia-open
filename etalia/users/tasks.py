# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.contrib.auth import get_user_model
from django.db.models import Count, F
from django.conf import settings
from pytz import timezone
from datetime import datetime, timedelta
from celery.canvas import chain
from mendeley.exception import MendeleyApiException

from config.celery import celery_app as app
from etalia.feeds.tasks import update_stream, update_trend, update_threadfeed
from etalia.popovers.tasks import init_popovers
from etalia.usersession.models import UserSession
from .models import UserLib, UserPeriodicEmail
from .constants import USERLIB_IDLE, USERLIB_NEED_REAUTH
from .utils import send_periodic_recommendation_email, next_weekday, \
    send_welcome_email

logger = logging.getLogger(__name__)

User = get_user_model()


@app.task(ignore_result=True)
def userlib_update_all():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .filter(lib__state=USERLIB_IDLE)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_lib.delay(user_pk)


@app.task(bind=True)
def update_lib(self, user_pk, **kwargs):
    """Async task for updating user library"""
    userlib = UserLib.objects.get(user_id=user_pk)
    try:
        userlib.update(**kwargs)
    except MendeleyApiException as exc:
        if not self.request.retries > 1:
            raise self.retry(exc=exc, countdown=60 * 10)
        else:
            # Remove session so that user need to re-authenticate
            userlib.set_state(USERLIB_NEED_REAUTH)
            UserSession.delete_user_sessions(self.user_id)
    return user_pk


@app.task()
def init_step(user_pk, step):
    """Set flag is_init to True"""
    user = User.objects.get(pk=user_pk)
    user.init_step = step
    user.save()

    return user_pk


def init_user(user_pk):
    """Task init user / Chain user library update, and feed initialization
    """
    task = chain(
        init_step.s(user_pk, 'LIB'),
        update_lib.s(),
        init_step.s('STR'),
        update_stream.s(),
        init_step.s('TRE'),
        update_trend.s(),
        init_step.s('THR'),
        update_threadfeed.s(),
        init_step.s('POP'),
        init_popovers.s(),
        init_step.s('IDL'),
    )

    task()
    return user_pk


@app.task(ignore_result=True)
def send_recommendation_emails_on_wed_11am():
    """Task that send batch of emails recommdation.

    Emails are sent to user on Wednesday 11am (in user timezone). Frequency varies
    depending on user settings.

    Task must be run once a week not to closed from the Monday 7am date +/- all
    timezones...
    """

    if settings.RECOMMENDATIONS_EMAILS_ON:
        us = User.objects.filter(settings__email_digest_frequency__gt=0)\
            .annotate(social_count=Count(F('social_auth')))\
            .exclude(social_count__lt=1)
        us = us.select_related('settings')
        us = us.select_related('userperiodicemail')

        # Get some date data
        utc = timezone('UTC')   # servers are in UTC
        nm = next_weekday(datetime.now().date(), 2)  # next wednesday

        for user in us:
            if hasattr(user, 'userperiodicemail') and user.userperiodicemail.last_sent_on:
                now = utc.localize(datetime.now())
                delta_since_last_email = (now - user.userperiodicemail.last_sent_on)
                user_period = user.settings.email_digest_frequency

                if timedelta(days=user_period) - delta_since_last_email > timedelta(days=5):
                    break

            # Do the math for the email to be send on monday 7am
            user_timezone = timezone(user.timezone)
            loc_dt = user_timezone.localize(datetime(nm.year, nm.month, nm.day, 11, 0, 0))
            utc_dt = loc_dt.astimezone(utc)
            # fire task
            async_send_periodic_email_at_eta.apply_async(args=[user.id, ],
                                                   eta=utc_dt)


@app.task(ignore_result=True)
def async_send_periodic_email_at_eta(user_id):
    """Send recommendations email"""
    send_periodic_recommendation_email(user_id)
    # Record email has been sent
    upe, _ = UserPeriodicEmail.objects.get_or_create(user_id=user_id)
    utc = timezone('UTC')
    upe.last_sent_on = utc.localize(datetime.now())
    upe.save()


@app.task(ignore_result=True)
def async_send_welcome_email(user_id):
    """Send welcoming email"""
    send_welcome_email(user_id)
