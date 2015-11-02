# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

FEED_STATUS_CHOICES = (
    ('NON', 'Uninitialized'),
    ('IDL', 'Idle'),
    ('ING', 'Syncing'),
)

STREAM_METHODS = (
    (0, 'Max'),
    (1, 'Average'),
    (2, 'Average Threshold binary'),
    (3, 'Average weighted journal'),
    (4, 'Average weighted journal-date'),
)

STREAM_METHODS_MAP = (
    (0, 'SimpleMax'),
    (1, 'SimpleAverage'),
    (2, 'ThresholdAverage'),
    (3, 'WeightedJournalAverage'),
    (4, 'WeightedJournalCreatedDateAverage'),
)

TREND_METHODS = (
    (1, 'Max'),
)


