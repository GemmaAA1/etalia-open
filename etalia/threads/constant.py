# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


THREAD_QUESTION = 1
THREAD_PAPER = 2
THREAD_TYPES = (
    (THREAD_QUESTION, 'Question'),
    (THREAD_PAPER, 'Paper'),
)


THREADFEED_STATUS_CHOICES = (
    ('NON', 'Uninitialized'),
    ('IDL', 'Idle'),
    ('ING', 'Syncing'),
)


THREAD_TIME_LAPSE_CHOICES = (
    (7, 'Week'),
    (30, 'Month'),
    (60, 'Two Months'),
    (180, 'Six Months'),
    (365, 'Year'),
)

# Invite
INVITE_PENDING = 1
INVITE_ACCEPTED = 2
INVITE_DECLINED = 3
INVITE_STATUSES = (
    (INVITE_PENDING, 'Pending'),
    (INVITE_ACCEPTED, 'Accepted'),
    (INVITE_DECLINED, 'Declined'),
)
