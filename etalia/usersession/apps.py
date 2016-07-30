# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.apps import AppConfig


class UserSessionConfig(AppConfig):
    name = 'etalia.usersession'
    verbose_name = "User Session"

    def ready(self):
        import etalia.usersession.signals  #noqa
