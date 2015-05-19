import os
from django.conf import settings

PUBLISHER_OPTIONS = [
    {'source': 'all',
     'file_path': os.path.join(str(settings.STATIC_ROOT),
                 'populate/publishers/publisher_list.csv')},
]

JOURNAL_OPTIONS = [
    {'source': 'thomson',
     'file_path': os.path.join(str(settings.STATIC_ROOT),
                 'populate/journals/20150510_thomsonreuters_cleaned.csv')},
    {'source': 'pubmed',
     'file_path': os.path.join(str(settings.STATIC_ROOT),
                 'populate/journals/20150510_medline_cleaned.csv')},
    {'source': 'arxiv',
     'file_path': os.path.join(str(settings.STATIC_ROOT),
                 'populate/journals/20150510_arxiv_cleaned.csv')},
]

# Consumer population is either defined through a file (e.g. pubmed) or by
# their publisher.
CONSUMER_OPTIONS = [
    {'source': 'pubmed',
     'file_path': os.path.join(str(settings.STATIC_ROOT),
                 'populate/journals/20150510_medline_cleaned.csv')},
]