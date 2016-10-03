# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

CONSUMER_THREAD_TYPE = (
    ('PPR', 'PubPeer'),
)

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
THREAD_INVITE_PENDING = 1
THREAD_INVITE_ACCEPTED = 2
THREAD_INVITE_DECLINED = 3
THREAD_INVITE_STATUSES = (
    (THREAD_INVITE_PENDING, 'Pending'),
    (THREAD_INVITE_ACCEPTED, 'Accepted'),
    (THREAD_INVITE_DECLINED, 'Declined'),
)

# PRIVACY
THREAD_PUBLIC = 1
THREAD_PRIVATE = 2
THREAD_PRIVACIES = (
    (THREAD_PUBLIC, 'Public'),
    (THREAD_PRIVATE, 'Private'),
)

# WATCH
THREAD_PINNED = 1
THREAD_BANNED = 2
THREAD_WATCH = (
    (THREAD_PINNED, 'Pinned'),
    (THREAD_BANNED, 'Banned'),
)

# PARTICIPATE
THREAD_JOINED = 1
THREAD_LEFT = 2
THREAD_PARTICIPATE = (
    (THREAD_JOINED, 'Joined'),
    (THREAD_LEFT, 'Left'),
)

# THIRD PARTY
# NB: Second member of tuple must match model name
THIRD_PARTY_PUBPEER = 1
THIRD_PARTY_TYPES = (
    (THIRD_PARTY_PUBPEER, 'PubPeer'),
)