# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import yaml
import copy
import random

from paperstream.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
from paperstream.nlp.models import Model, MostSimilar
from paperstream.library.models import Paper
from paperstream.altmetric.tasks import update_altmetric
from paperstream.altmetric.models import AltmetricModel

# Inputs
MODEL_FILE = './routines/models.yaml'
MODEL_NAME = 'test'

# POPULATE DATABASE WITH SOME DATA

# fetching some recent papers from the consumers
pubmed_run('pubmed_all')
arxiv_run('arxiv_all')
elsevier_run('elsevier_all')

# Build NLP model and activate
with open(MODEL_FILE) as stream:
    MODELS = yaml.load(stream)
model_params = [m for m in MODELS if m['name'] == MODEL_NAME][0]

model = Model.objects.create(**model_params)
papers = Paper.objects.filter(is_trusted=True,
                              title__regex=r'^.{5}.*',
                              abstract__regex=r'^.{10}.*')
model.dump(papers)
model.build_vocab_and_train()
model.activate()

# Populate paper with NLP signatures
model.save_journal_vec_from_bulk()
model.save_paper_vec_from_bulk()

# Build MostSimilar for paper matching
ms, _ = MostSimilar.objects.get_or_create(model=model)
ms.full_update()
ms.activate()

# Fetch some Altmetric data and fake others (to save time)
alt_objs = []
for paper in papers[:10]:
    update_altmetric(paper.pk)
    obj = copy.copy(paper.altmetric.__dict__)
    obj.pop('paper_id')
    obj.pop('id')
    obj.pop('_state')
    obj.pop('_paper_cache')
    alt_objs.append(obj)
# and randomly assigned to other to save time
for paper in papers[10:]:
    alt_obj = alt_objs[random.randint(0, 9)]
    altmetric, _ = AltmetricModel.objects.get_or_create(paper_id=paper.pk,
                                                        **alt_obj)


