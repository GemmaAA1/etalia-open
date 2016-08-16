# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib.auth.signals import user_logged_in
from django.dispatch import receiver
from easy_timezones.signals import detected_timezone
from .models import UserSession
from django.contrib.auth import get_user_model

User = get_user_model()


@receiver(user_logged_in)
def user_logged_in_handler(sender, request, user, **kwargs):
    if request.session.session_key:
        UserSession.objects.get_or_create(
            user=user,
            session_id=request.session.session_key
        )


@receiver(detected_timezone, sender=User)
def process_timezone(sender, instance, timezone, **kwargs):
    if instance.timezone != timezone:
        instance.timezone = timezone
        instance.save()
