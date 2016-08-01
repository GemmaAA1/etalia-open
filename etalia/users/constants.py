# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# Define constants used in etalia
USER_INDIVIDUAL = 1
USER_THIRD_PARTY = 2
USER_TYPES = (
    (USER_INDIVIDUAL, 'Individual'),
    (USER_THIRD_PARTY, 'Third-Party')
)


# Type of providers (defined in users app)
PROVIDER_TYPE = (
    ('ZOT', 'Zotero'),
    ('MEN', 'Mendeley'),
    ('', 'Unknown')
)

MENDELEY_PT = (
    ('journal',                 'JOU'),
    ('book',                    'BOO'),
    ('book_section',            'BOS'),
    ('conference_proceedings',  'PRO'),
    ('working_paper',           'DRA'),
    ('patent',                  'PAT'),
    ('thesis',                  'THE'),
    ('UNKNOWN',                 ''),
)

ZOTERO_PT = (
    ('journalArticle',  'JOU'),
    ('book',            'BOO'),
    ('bookSection',     'BOS'),
    ('conferencePaper', 'PRO'),
    ('thesis',          'THE'),
    ('patent',          'PAT'),
    ('letter',          'LET'),
    ('UNKNOWN',         ''),
)

INIT_STEPS = (
    ('NON', 'uninitialized'),
    ('LIB', 'Syncing library'),
    ('STR', 'Syncing feed Papers'),
    ('TRE', 'Syncing feed Trend'),
    ('THR', 'Syncing feed Thread'),
    ('POP', 'Syncing popovers'),
    ('IDL', 'done'),
)

# Relationship
RELATIONSHIP_FOLLOWING = 1
RELATIONSHIP_BLOCKED = 2
RELATIONSHIP_STATUSES = (
    (RELATIONSHIP_FOLLOWING, 'Following'),
    (RELATIONSHIP_BLOCKED, 'Blocked'),
)

# UserLib
USERLIB_UNINITIALIZED = 1
USERLIB_IDLE = 2
USERLIB_SYNCING = 3
USERLIB_NEED_REAUTH = 4
USERLIB_STATE_CHOICES = (
    (USERLIB_UNINITIALIZED, 'Uninitialized'),
    (USERLIB_IDLE, 'Idle'),
    (USERLIB_SYNCING, 'Syncing'),
    (USERLIB_NEED_REAUTH, 'Need New OAuth')
)