# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import redirect
from django.conf import settings
from social.pipeline.partial import partial

from paperstream.core.utils import get_celery_worker_status
from .models import Affiliation

from .tasks import update_lib as async_update_lib
from .tasks import init_user as async_init_user


@partial
def require_primary(strategy, details, user=None, is_new=False, *args, **kwargs):
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
def require_affiliation(strategy, details, request=None, user=None, *args, **kwargs):
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
def init_user(social, user, *args, **kwargs):
    async_init_user(user.pk, social.provider)
    user.lib.set_state('ING')
    user.feeds.first().set_state('ING')
    return {}


