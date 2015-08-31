# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

from .models import User, UserSettings, UserLib

admin.site.register(User)
admin.site.register(UserLib)
admin.site.register(UserSettings)

