# Define constants used in paperstream


# Type of consumer (defined in consumers app)
CONSUMER_TYPE = (
    ('PUMD', 'PubMed'),
    ('ELSV', 'Elsevier'),
    ('ARXI', 'Arxiv'),
    ('IEEE', 'IEEE'),
)

# Pubmed publication matching
PUBMED_PT = (
    ('Journal Article', 'JOUR'),
    ('Letter',          'LETT'),
    ('Editorial',       'EDIT'),
    ('News',            'NEWS'),
    ('Congresses',      'PROC'),
    ('Review',          'REVI'),
    ('Unknown',         ''),
)