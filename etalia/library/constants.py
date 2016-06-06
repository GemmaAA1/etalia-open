# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from etalia.consumers.constants import CONSUMER_TYPE
from etalia.users.constants import PROVIDER_TYPE

# Define constants used in models

SOURCE_TYPE = tuple(set(CONSUMER_TYPE + PROVIDER_TYPE))

PAPER_TYPE = (
    ('JOU', 'Journal Article'),
    ('LET', 'Letter'),
    ('EDI', 'Editorial'),
    ('NEW', 'News'),
    ('PRO', 'Proceedings'),
    ('REV', 'Review'),
    ('PRE', 'e-Print'),  # e.g Arxiv
    ('DRA', 'Draft'),
    ('BOO', 'Book'),
    ('BOS', 'Book section'),
    ('PAT', 'Patent'),
    ('THE', 'Thesis'),
    ('', 'Unknown')
)

PUBLISH_STATUS = (
    ('ppublish', 'Paper Print'),
    ('epublish', 'e-Print only'),
    ('aheadofprint', 'e-Print ahead'),
    ('preprint', 'pre-Print'),
    ('', 'Unknown')
)

# Used in Journal
PUBLISH_PERIODS = (
    ('ANN', 'Annual'),
    ('SEM', 'Semi-annual'),
    ('TRI', 'Tri-annual'),
    ('QUA', 'Quarterly'),
    ('MON', 'Monthly'),
    ('BIM', 'Bi-monthly'),
    ('IRR', 'Irregular'),
    ('', 'Unknown'))

# Language code used in Journal, Paper
# NB: languages code are from pubmed
LANGUAGES = (
    ('ENG', 'English'),
    ('AFR', 'Afrikaans'),
    ('ALB', 'Albanian'),
    ('AMH', 'Amharic'),
    ('ARA', 'Arabic'),
    ('ARM', 'Armenian'),
    ('AZE', 'Azerbaijani'),
    ('BEN', 'Bengali'),
    ('BOS', 'Bosnian'),
    ('BUL', 'Bulgarian'),
    ('CAT', 'Catalan'),
    ('CHI', 'Chinese'),
    ('CZE', 'Czech'),
    ('DAN', 'Danish'),
    ('DUT', 'Dutch'),
    ('ENG', 'English'),
    ('EPO', 'Esperanto'),
    ('EST', 'Estonian'),
    ('FIN', 'Finnish'),
    ('FRE', 'French'),
    ('GEO', 'Georgian'),
    ('GER', 'German'),
    ('GLA', 'Scottish Gaelic'),
    ('GRE', 'Greek, Modern'),
    ('HEB', 'Hebrew'),
    ('HIN', 'Hindi'),
    ('HRV', 'Croatian'),
    ('HUN', 'Hungarian'),
    ('ICE', 'Icelandic'),
    ('IND', 'Indonesian'),
    ('ITA', 'Italian'),
    ('JPN', 'Japanese'),
    ('KIN', 'Kinyarwanda'),
    ('KOR', 'Korean'),
    ('LAT', 'Latin'),
    ('LAV', 'Latvian'),
    ('LIT', 'Lithuanian'),
    ('MAC', 'Macedonian'),
    ('MAL', 'Malayalam'),
    ('MAO', 'Maori'),
    ('MAY', 'Malay'),
    ('MUL', 'Multiple languages'),
    ('NOR', 'Norwegian'),
    ('PER', 'Persian, Iranian'),
    ('POL', 'Polish'),
    ('POR', 'Portuguese'),
    ('RUM', 'Romanian, Rumanian, Moldovan'),
    ('RUS', 'Russian'),
    ('SAN', 'Sanskrit'),
    ('SLO', 'Slovak'),
    ('SLV', 'Slovenian'),
    ('SPA', 'Spanish'),
    ('SRP', 'Serbian'),
    ('SWE', 'Swedish'),
    ('THA', 'Thai'),
    ('TUR', 'Turkish'),
    ('UKR', 'Ukrainian'),
    ('UND', 'Undetermined'),
    ('URD', 'Urdu'),
    ('VIE', 'Vietnamese'),
    ('WEL', 'Welsh')
)

# Language code 2 pap code from language detection
LANGUAGES_DETECT = (
    ('AF', 'AFR'),
    ('SQ', 'ALB'),
    ('AM', 'AMH'),
    ('AR', 'ARA'),
    ('HY', 'ARM'),
    ('AZ', 'AZE'),
    ('BN', 'BEN'),
    ('BS', 'BOS'),
    ('BG', 'BUL'),
    ('CA', 'CAT'),
    ('ZH', 'CHI'),
    ('CS', 'CZE'),
    ('DA', 'DAN'),
    ('NL', 'DUT'),
    ('EN', 'ENG'),
    ('EO', 'EPO'),
    ('ET', 'EST'),
    ('FI', 'FIN'),
    ('FR', 'FRE'),
    ('KA', 'GEO'),
    ('DE', 'GER'),
    ('GD', 'GLA'),
    ('EL', 'GRE'),
    ('IW', 'HEB'),
    ('HI', 'HIN'),
    ('HR', 'HRV'),
    ('HU', 'HUN'),
    ('IS', 'ICE'),
    ('ID', 'IND'),
    ('IT', 'ITA'),
    ('JA', 'JPN'),
    ('RW', 'KIN'),
    ('KO', 'KOR'),
    ('LA', 'LAT'),
    ('LV', 'LAV'),
    ('LT', 'LIT'),
    ('MK', 'MAC'),
    ('ML', 'MAL'),
    ('MI', 'MAO'),
    ('MS', 'MAY'),
    ('NO', 'NOR'),
    ('FA', 'PER'),
    ('PL', 'POL'),
    ('PT', 'POR'),
    ('RO', 'RUM'),
    ('RU', 'RUS'),
    ('HI', 'SAN'),
    ('SK', 'SLO'),
    ('SL', 'SLV'),
    ('ES', 'SPA'),
    ('SR', 'SRP'),
    ('SV', 'SWE'),
    ('TH', 'THA'),
    ('TR', 'TUR'),
    ('UK', 'UKR'),
    ('UR', 'URD'),
    ('VI', 'VIE'),
    ('CY', 'WEL'),
    (''  , 'UND'),
)


PAPER_TIME_LAPSE_CHOICES = (
    (7, 'Week'),
    (30, 'Month'),
    (60, 'Two Months'),
    (180, 'Six Months'),
    (365, 'Year'),
)


# WATCH
PAPER_PINNED = 1
PAPER_BANNED = 2
PAPER_WATCH = (
    (PAPER_PINNED, 'Pinned'),
    (PAPER_BANNED, 'Banned'),
)

# ADD
PAPER_ADDED = 1
PAPER_TRASHED = 2
PAPER_STORE = (
    (PAPER_ADDED, 'Pinned'),
    (PAPER_TRASHED, 'Trashed'),
)
