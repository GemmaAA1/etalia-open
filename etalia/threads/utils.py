# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.contrib.auth import get_user_model
from django.conf import settings
from etalia.core.emails import Email
from .models import Thread

User = get_user_model()


def send_invite_thread_email(from_id, to_id, thread_id):
    """Send etalia welcome email"""
    from_user = User.objects.get(id=from_id)
    to_user = User.objects.get(id=to_id)
    thread = Thread.objects.get(id=thread_id)

    # Instantiate
    email = Email(
        template=settings.INVITE_THREAD_EMAIL_TEMPLATE,
        cids={},
        tags=['invite-thread'],
        metadata={'from_user': from_user.id,
                  'to_user': to_user.id},
        subject='{0} invites you to join a thread'.format(from_user.get_full_name()),
        from_email='etalia@etalia.org',
        to=[to_user.email],
        reply_to=['contact@etalia.org'],
        extra_ctx={'thread': thread,
                   'from_user': from_user.get_full_name(),
                   'to_user': to_user.get_full_name()})

    # send
    email.send()
