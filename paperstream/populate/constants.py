# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
from django.conf import settings

PUBLISHER_OPTIONS = [
    {'source': 'all',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/publishers/publisher_list.csv')},
]

JOURNAL_OPTIONS = [
    {'source': 'thomson',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_thomsonreuters_cleaned.csv')},
    {'source': 'pubmed',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_medline_cleaned.csv')},
    {'source': 'arxiv',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_arxiv_cleaned.csv')},
    {'source': 'thomson_local',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_thomsonreuters_cleaned_light.csv')},
    {'source': 'pubmed_local',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_medline_cleaned_light.csv')},
    {'source': 'arxiv_local',
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_arxiv_cleaned_light.csv')},
]

# Consumer population is either defined through a file (e.g. pubmed) or by
# their publisher.
CONSUMER_OPTIONS = [
    {'source': 'pubmed',
     'local': False,
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_medline_cleaned.csv')},
    {'source': 'pubmed',
     'local': True,
     'file_path': os.path.join(str(settings.STATICFILES_DIRS[0]),
                 'populate/journals/20150510_medline_cleaned_light.csv')},
]