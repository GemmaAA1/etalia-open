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