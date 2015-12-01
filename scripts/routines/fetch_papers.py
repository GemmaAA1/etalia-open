# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import yaml

from paperstream.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
from paperstream.nlp.models import Model, MostSimilar
from paperstream.library.models import Paper

# Inputs
MODEL_FILE = './routines/models.yaml'
MODEL_NAME = 'test'

# fetching some papers
pubmed_run('pubmed_all')
arxiv_run('arxiv_all')
elsevier_run('elsevier_all')

# building NLP model
with open(MODEL_FILE) as stream:
    MODELS = yaml.load(stream)
model_params = [model for model in MODELS if model['name'] == MODEL_NAME]

model = Model.objects.create(**model_params)
papers = Paper.objects.filter(is_trusted=True,
                              title__regex=r'^.{5}.*',
                              abstract__regex=r'^.{10}.*')
model.dump(papers)
model.build_vocab_and_train()

# Propagate to MostSimilar, journalvector, papervector
model.propagate()
