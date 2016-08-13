# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models
from django.conf import settings
from django.contrib.sessions.models import Session

from etalia.core.models import TimeStampedModel


class UserSession(TimeStampedModel):
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    session = models.ForeignKey(Session)

    @classmethod
    def delete_user_sessions(cls, user_id):
        user_sessions = cls.objects.filter(user_id=user_id)
        for user_session in user_sessions:
            user_session.session.flush()
