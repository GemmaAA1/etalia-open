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

LSHF_STATES = (
    ('NON', 'None'),
    ('BUS', 'Busy'),
    ('IDL', 'Idle'),
)
