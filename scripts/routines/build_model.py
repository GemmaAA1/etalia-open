# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import yaml
from paperstream.nlp.models import Model, LSH
from paperstream.library.models import Paper

path_to_model_file = 'models.yaml'
with open(path_to_model_file) as stream:
    MODELS = yaml.load(stream)


def build(model_name, papers=None):

    check_unique_name(MODELS)
    model_args = [model for model in MODELS if model['name'] == model_name][0]

    # Initiate model
    model = Model.objects.create(**model_args)
    # dump papers data
    if not papers:
        papers = Paper.objects.all()
    model.dump(papers)
        model.build_vocab_and_train()
    # Propagate to LSH, journalvector, papervector
    model.propagate()


def check_unique_name(models):
    names = [(i, models['name']) for i, model in enumerate(models)]
    if not len(names) == len(set(names)):
        raise AssertionError('MODELS has duplicated names')