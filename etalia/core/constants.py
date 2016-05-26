# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf import settings

EMAIL_DIGEST_FREQUENCY_CHOICES = (
    (7, 'Weekly'),
    (15, 'Bi-Weekly'),
    (30, 'Monthly'),
    (30, 'Bi-Monthly'),
    (-1, 'Never'),
)