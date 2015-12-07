# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

FEED_STATUS_CHOICES = (
    ('NON', 'Uninitialized'),
    ('IDL', 'Idle'),
    ('ING', 'Syncing'),
)

STREAM_METHODS = (
    (0, 'Occurrence Count'),
    (1, 'Max'),
    (2, 'Average'),
    (3, 'Average Threshold binary'),
    (4, 'Average weighted journal'),
    (5, 'Average weighted journal-date'),
)

STREAM_METHODS_MAP = (
    (0, 'OccurrenceCount'),
    (1, 'SimpleMax'),
    (2, 'SimpleAverage'),
    (3, 'ThresholdAverage'),
    (4, 'WeightedJournalAverage'),
    (5, 'WeightedJournalCreatedDateAverage'),
)

TREND_METHODS = (
    (0, 'Method #1'),
)


