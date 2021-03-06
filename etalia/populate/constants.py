# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
from django.conf import settings

PUBLISHER_OPTIONS = [
    {'source': 'all',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/publishers/publisher_list.json')},
]

JOURNAL_OPTIONS = [
    {'source': 'thomson',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_thomsonreuters_cleaned.json')},
    {'source': 'pubmed',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_medline_cleaned.json')},
    {'source': 'pubmed2',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20160225_Entrez_parse.json')},
    {'source': 'arxiv',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_arxiv_cleaned.json')},
    {'source': 'thomson_local',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_thomsonreuters_cleaned_light.json')},
    {'source': 'pubmed_local',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_medline_cleaned_light.json')},
    {'source': 'arxiv_local',
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_arxiv_cleaned_light.json')},
]

# Consumer population is either defined through a file (e.g. pubmed) or by
# their publisher.
CONSUMER_OPTIONS = [
    {'source': 'pubmed',
     'local': False,
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20160225_Entrez_parse.json')},
    {'source': 'pubmed',
     'local': True,
     'file_path': os.path.join(str(settings.ROOT_DIR), 'data',
                 'populate/journals/20150510_medline_cleaned_light.json')},
]