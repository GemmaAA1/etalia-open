# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models


class Press(models.Model):

    title = models.CharField(max_length=256)

    body = models.TextField()

    location = models.CharField(max_length=256)

    date = models.DateField(auto_now_add=True)

    canonical_url = models.URLField(null=True, blank=True, default='')