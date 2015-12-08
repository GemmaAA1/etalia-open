# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf import settings

NLP_TIME_LAPSE_CHOICES = (
    (7, 'Week'),
    (30, 'Month'),
    (60, 'Two Months'),
    (180, 'Six Months'),
    (365, 'Year'),
)

NLP_NARROWNESS_CHOICES = (
    (-1, 'Narrower'),
    (0, 'Normal'),
    (1, 'Broader'),
)

EMAIL_DIGEST_FREQUENCY_CHOICES = (
    (7, 'Weekly'),
    (15, 'Bi-Weekly'),
    (30, 'Monthly'),
    (30, 'Bi-Monthly'),
    (-1, 'Never'),
)