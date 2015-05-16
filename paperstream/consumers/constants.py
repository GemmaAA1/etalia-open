# Define constants used in paperstream


# Type of consumer (defined in consumers app)
CONSUMER_TYPE = (
    ('PUBM', 'PubMed'),
    ('ELSV', 'Elsevier'),
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
    ('UBKNOWN',         ''),
)