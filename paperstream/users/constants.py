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
     'Below the stream of the most relevant scientific articles we found for ' \
     'you. You can filter them using the right-side panel or using the shortcut icons in the top navigation '
     'bar. Pins (<span class="eai eai-pin"></span>) are used to quickly mark the articles'
     'you are interested in. Bans (<span class="eai eai-remove"></span>) will remove the article '
     'from your stream. Using <span class="eai eai-pin"></span> and '
     '<span class="eai eai-remove"></span> allows etalia to better understand your interests '
     'and improve your stream content.'),
    (2,
     'trend',
     'Below a list of publications that are currently receiving the most ' \
     'attention. You can filter them using the right-side panel or using the shortcut icons in the top navigation '
     'bar. Pins (<span class="eai eai-pin"></span>) are used to quickly mark the articles'
     'you are interested in. Bans (<span class="eai eai-remove"></span>) will remove the article '
     'from your stream. Using <span class="eai eai-pin"></span> and '
     '<span class="eai eai-remove"></span> allows etalia to better understand your interests '
     'and improve your trend content.')
)
