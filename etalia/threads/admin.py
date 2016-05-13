# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

from .models import Thread, ThreadUserInvite, ThreadUser, ThreadPost, \
    ThreadComment

admin.site.register(Thread)
admin.site.register(ThreadUser)
admin.site.register(ThreadUserInvite)
admin.site.register(ThreadPost)
admin.site.register(ThreadComment)
