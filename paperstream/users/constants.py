# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# Define constants used in paperstream

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
    ('LIB', 'library'),
    ('STR', 'personalized stream'),
    ('TRE', 'trends'),
    ('IDL', 'done'),
)

# MESSAGES
# User linked persistant message (i.e. that user needs to close)
# in the form (id, extra_tag (str), message (str))
PERSISTANT_USER_MESSAGES = (
    (1,
     'stream',
     'You find here the most relevant scientific publications we found for ' \
     'you based on the publications that you stored in your reference manager.'),
    (2,
     'trend',
     'You find here a list of publications that are currently receiving the most ' \
     'attention. A new way to broaden your interests and get aware of potential breakthroughs before other')
)
