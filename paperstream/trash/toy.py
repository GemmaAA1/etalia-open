# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from paperstream.nlp.models import Model
model_names = ['dbow-64-mc-5', 'dbow-128-with-words', 'dm-128']

for model_name in model_names:
    model = Model.objects.load(name=model_name)
    model.propagate()

model_name = 'dbow-128'
import yaml

path_to_model_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/scripts/routines/models.yaml'
with open(path_to_model_file) as stream:
    MODELS = yaml.load(stream)
model_args = [m for m in MODELS if m['name'] == model_name][0]
model = Model.objects.create(**model_args)
model.build_vocab_and_train()
# Propagate to LSH, journalvector, papervector
model.propagate()



from django.db import transaction
import time
from paperstream.nlp.models import LSH

lsh = LSH.objects.load(model__name='dbow-64-mc-5', time_lapse=-1)
t0 = time.time()
pks = lsh.lsh.pks[1000:1100]
lsh.update_neighbors(pks)
dt = time.time() - t0
print(dt)
