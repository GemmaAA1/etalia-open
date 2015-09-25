# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import yaml
from paperstream.nlp.models import Model
from paperstream.library.models import Paper

path_to_model_file = 'models.yaml'
with open(path_to_model_file) as stream:
    MODELS = yaml.load(stream)


def build(model_name, papers=None):

    model_args = [m for m in MODELS if m['name'] == model_name][0]

    # Initiate model
    model = Model.objects.create(**model_args)
    # dump papers data
    if not papers:
        papers = Paper.objects.filter(is_trusted=True,
                                      title__regex=r'^.{5}.*',
                                      abstract__regex=r'^.{10}.*')
    model.dump(papers)
    model.build_vocab_and_train()
    # Propagate to MostSimilar, journalvector, papervector
    model.propagate()
