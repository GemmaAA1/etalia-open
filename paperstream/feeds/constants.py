# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

FEED_STATUS_CHOICES = (
    ('NON', 'Uninitialized'),
    ('IDL', 'Idle'),
    ('ING', 'Syncing'),
)


FEED_SCORING_CHOICES = (
    (0, 'Max'),
    (1, 'Average'),
    (2, 'Average Threshold binary'),
    (3, 'Average weighted journal'),
    (4, 'Average weighted journal-date'),
)
