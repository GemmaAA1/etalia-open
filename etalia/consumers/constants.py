# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# Define constants used in etalia


# Type of consumer (defined in consumers app)
# WARNING: HUMAN READABLE MUST MATCH Publisher.name FOR CONSUMER POPULATION
# MANAGEMENT ROUTINE
CONSUMER_TYPE = (
    ('PUB', 'PubMed'),
    ('ELS', 'Elsevier'),
    ('ARX', 'Arxiv'),
    ('IEE', 'IEEE'),
)
