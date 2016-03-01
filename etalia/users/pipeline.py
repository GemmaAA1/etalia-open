# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import urllib.request
from urllib.parse import urlparse
from django.core.files import File
from django.shortcuts import redirect
from django.conf import settings
from social.pipeline.partial import partial
from avatar.models import Avatar

from etalia.core.utils import get_celery_worker_status
from messages_extends.models import Message
from messages_extends import constants as constants_messages
from .models import Affiliation

from .forms import UserAffiliationForm
from .tasks import update_lib as async_update_lib
from .tasks import init_user as async_init_user
from .constants import PERSISTANT_USER_MESSAGES

@partial
def require_primary(strategy, details, *args, user=None, **kwargs):
    """ Redirect to primary info form for user to check them
    """
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

@partial
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
    if photo_url:
        url_parsed = urlparse(photo_url)
        photo_filename = os.path.basename(url_parsed.path)
        # if photo not in [settings.AVATAR_DEFAULT_MENDELEY, settings.AVATAR_DEFAULT_ZOTERO]:
        local_filename, headers = urllib.request.urlretrieve(photo_url)
        f = open(local_filename, 'rb')
        avatar = Avatar(user=user, primary=True)
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


@partial
def update_user_lib(backend, social, user, *args, **kwargs):
    session = backend.get_session(social, user)
    if settings.DEBUG and get_celery_worker_status().get('ERROR'):
        backend.update_lib(user, session)
    else:  # celery is running -> go async
        async_update_lib.apply_async(args=[user.pk, social.provider],
                                     serializer='json')
    return {}


@partial
def init_messages(social, user, *args, **kwargs):

    for id, extra_tags, message in PERSISTANT_USER_MESSAGES:
        mess, new = Message.objects.get_or_create(
            user=user,
            level=constants_messages.INFO_PERSISTENT,
            extra_tags='id{id},{tags}'.format(id=id, tags=extra_tags))
        if new:
            mess.message = message
            mess.save()

    return {}

@partial
def init_user(social, user, *args, **kwargs):
    if user.lib.state == 'NON':  # user non-initialized yet
        async_init_user(user.pk, social.provider)
    return {}


