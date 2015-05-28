# Define constants used in paperstream


# Type of consumer (defined in consumers app)
# WARNING: HUMAN READABLE MUST MATCH Publisher.name FOR CONSUMER POPULATION
# MANAGEMENT ROUTINE
CONSUMER_TYPE = (
    ('PUB', 'PubMed'),
    ('ELS', 'Elsevier'),
    ('ARX', 'Arxiv'),
    ('IEE', 'IEEE'),
)

# Pubmed publication matching
PUBMED_PT = (
    ('JOURNAL ARTICLE', 'JOU'),
    ('LETTER',          'LET'),
    ('EDITORIAL',       'EDI'),
    ('NEWS',            'NEW'),
    ('CONGRESSES',      'PRO'),
    ('REVIEW',          'REV'),
    ('PATENTS',         'PAT'),
    ('UNKNOWN',         ''),
)