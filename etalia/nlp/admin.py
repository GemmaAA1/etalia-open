# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

# Register your models here.

from .models import Model, PaperEngine

admin.site.register(Model)
admin.site.register(PaperEngine)