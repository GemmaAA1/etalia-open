# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

FEED_STATUS_CHOICES = (
    ('NON', 'Uninitialized'),
    ('IDL', 'Idle'),
    ('ING', 'Syncing'),
)

STREAM_METHODS = (
    (0, 'Content Based Simple'),
    # (1, 'Max'),
    # (2, 'Average'),
    # (3, 'Average Threshold binary'),
    # (4, 'Average weighted journal'),
    # (5, 'Average weighted journal-date'),
    # (6, 'Occurrence Count'),
)

STREAM_METHODS_MAP = (
    (0, 'ContentBasedSimple'),
    # (1, 'SimpleMax'),
    # (2, 'SimpleAverage'),
    # (3, 'ThresholdAverage'),
    # (4, 'WeightedJournalAverage'),
    # (5, 'WeightedJournalCreatedDateAverage'),
    # (6, 'OccurrenceCount'),
)

TREND_METHODS = (
    (0, 'Method #1'),
)


