# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

MODEL_STATES = (
    ('UNT', 'Untrained'),
    ('VOC', 'Building Vocabulary'),
    ('TRA', 'Training'),
    ('SAV', 'Saving'),
    ('LOA', 'Loading'),
    ('POP', 'Populating'),
    ('IDL', 'Idle'),
    ('USE', 'Usable'),
)

FIELDS_FOR_MODEL = (
    ('title', 'Title'),
    ('abstract', 'Abstract'),
)

NLP_JOURNAL_RATIO_CHOICES = (
    ('0.0', '0 %'),
    ('0.05', '5 %'),
    ('0.10', '10 %'),
    ('0.15', '15 %'),
    ('0.20', '20 %'),
    ('0.25', '25 %'),
    ('0.30', '30 %'),
)
