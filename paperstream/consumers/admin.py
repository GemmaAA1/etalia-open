# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

# Register your models here.

from .models import ConsumerElsevier, ConsumerArxiv, ConsumerPubmed

admin.register(ConsumerArxiv)
admin.register(ConsumerElsevier)
admin.register(ConsumerPubmed)