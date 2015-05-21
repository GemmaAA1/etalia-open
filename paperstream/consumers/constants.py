# Define constants used in paperstream


# Type of consumer (defined in consumers app)
# WARNING: HUMAN READABLE MUST MATCH Publisher.name FOR CONSUMER POPULATION
# MANAGEMENT ROUTINE
CONSUMER_TYPE = (
    ('PUBM', 'PubMed'),
    ('ELSE', 'Elsevier'),
    ('ARXI', 'Arxiv'),
    ('IEEE', 'IEEE'),
)

# Pubmed publication matching
PUBMED_PT = (
    ('JOURNAL ARTICLE', 'JOUR'),
    ('LETTER',          'LETT'),
    ('EDITORIAL',       'EDIT'),
    ('NEWS',            'NEWS'),
    ('CONGRESSES',      'PROC'),
    ('REVIEW',          'REVI'),
    ('UNKNOWN',         ''),
)