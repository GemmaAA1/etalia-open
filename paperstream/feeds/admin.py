# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

# Register your models here.

from .models import UserFeed


class PaperAdmin(admin.ModelAdmin):
    search_fields = ('user', 'name')
    list_display = ('user', 'name', 'state')

admin.register(UserFeed)
