# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib.auth.signals import user_logged_in
from .models import UserSession
from django.dispatch import receiver


@receiver(user_logged_in)
def user_logged_in_handler(sender, request, user, **kwargs):
    if request.session.session_key:
        UserSession.objects.get_or_create(
            user=user,
            session_id=request.session.session_key
        )
