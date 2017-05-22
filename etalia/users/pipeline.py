# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import urllib.request
from urllib.parse import urlparse
from django.core.files import File
from django.shortcuts import redirect
from django.core.mail import send_mail
from django.conf import settings
from django.contrib.staticfiles.templatetags.staticfiles import static
from social.pipeline.partial import partial
from avatar.models import Avatar
from etalia.usersession.models import UserSession
from .models import Affiliation
from .constants import USERLIB_UNINITIALIZED

from .forms import UserAffiliationForm
from .tasks import init_user as async_init_user, async_send_welcome_email


@partial
def require_primary(strategy, details, *args, user=None, **kwargs):
    """ Redirect to primary info form for user to check them"""
    if user and user.email and user.first_name and user.last_name:
        return
    else:
        if strategy.session.get('basic_info', None):
            basic_info = strategy.session.get('basic_info')
            details['first_name'] = basic_info.get('first_name')
            details['last_name'] = basic_info.get('last_name')
            details['email'] = basic_info.get('email')
            return
        return redirect('user:require-basic-info')


def create_details(strategy, details, *args, user=None, **kwargs):
    # link affiliation
    affiliation_kwargs = details.get('tmp_affiliation')
    if affiliation_kwargs:
        try:
            affiliation = Affiliation.objects.get(**affiliation_kwargs)
            user.affiliation = affiliation
        except Affiliation.DoesNotExist:
            form = UserAffiliationForm(affiliation_kwargs)
            if form.is_valid:
                affiliation = form.save()
                user.affiliation = affiliation
    # link avatar
    photo_url = details.get('photo')
    if not photo_url:
        photo_url = settings.AVATAR_DEFAULT_URL
    url_parsed = urlparse(photo_url)
    photo_filename = os.path.basename(url_parsed.path)
    local_filename, headers = urllib.request.urlretrieve(photo_url)
    avatar = Avatar(user=user, primary=True)
    with open(local_filename, 'rb') as f:
        avatar.avatar.save(photo_filename, File(f))
    avatar.save()

    # link title
    user.title = details.get('title', '')
    # link position
    user.position = details.get('position', '')

    user.save()

    return {}


@partial
def require_affiliation(strategy, details, *args, user=None, **kwargs):
    """Redirect to affiliation info form for user to check them"""
    if getattr(user, 'affiliation'):
        return
    else:
        if strategy.session.get('affiliation_pk', None):
            affiliation_pk = strategy.session.get('affiliation_pk')
            try:
                affiliation = Affiliation.objects.get(pk=affiliation_pk)
                user.affiliation = affiliation
                user.save()
                return
            except Affiliation.DoesNotExist:
                pass

        return redirect('user:require-affiliation')


def update_usersession(strategy, details, *args, **kwargs):
    user = kwargs.get('user')
    UserSession.objects.get_or_create(
        user=user,
        session_id=strategy.request.session.session_key)
    return {}


def send_email_of_new_signup(strategy, details, *args, **kwargs):
    emails = [u[1] for u in settings.ADMINS]
    user = kwargs.get('user')
    if not user.last_login:
        try:
            send_mail('New Signup ({0})'.format(user.email),
                       '{0} just signed-up'.format(user.email),
                       'etalia@etalia.io',
                        emails,
                       fail_silently=False)
        # TODO specify exceptions
        except:
            pass


def send_welcome_email_at_signup(strategy, details, *args, **kwargs):
    user = kwargs.get('user')
    if not user.last_login:
        try:
            async_send_welcome_email.delay(user.id)
        # TODO specify exceptions
        except:
            pass


def init_user(social, user, *args, **kwargs):
    if user.lib.state == USERLIB_UNINITIALIZED:  # user non-initialized yet
        async_init_user(user.pk)
    return {}




