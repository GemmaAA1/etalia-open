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
     'You find here a stream of the most relevant scientific articles we found for ' \
     'you based on the publications that you stored in your reference manager. You can'
     'filter them using the right-side panel or using the shortcut icons in the top navigation'
     'bar. Pins (<span class="eai eai-pin"></span>) are used to quickly mark the articles'
     'you are interested in before adding them to your library. Bans ('
     '<span class="eai eai-ban"></span>) will remove the article'
     'from your stream. Using <span class="eai eai-pin"></span> and '
     '<span class="eai eai-ban"></span> allow etalia to better understanding your interests'
     'and to provide you with ever improving content.'),
    (2,
     'trend',
     'You find here a list of publications that are currently receiving the most ' \
     'attention. A new way to broaden your interests and get aware of potential '
     'breakthroughs before your peers. You can'
     'filter them using the right-side panel or using the shortcut icons in the top navigation'
     'bar. Pins (<span class="eai eai-ban"></span>) are used to quickly mark the articles'
     'you are interested in before adding them to your library. Bans ('
     '<span class="eai eai-ban"></span>) will remove the article'
     'from your stream. Using <span class="eai eai-pin"></span> and '
     '<span class="eai eai-ban"></span> allow etalia to better understanding your interests'
     'and to provide you with ever improving content.')
)
