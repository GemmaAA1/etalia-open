# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import csv
import os
import numpy as np
from django.contrib.auth import get_user_model
from paperstream.library.models import Paper
from paperstream.nlp.models import PaperVectors

email = 'norbert.schuff@gmail.com'
email = 'nicolas.pannetier@gmail.com'
out_path = '~/.'

User = get_user_model()

user = User.objects.get(email=email)

library = user.lib.papers.all()

# Store vector user library
user_title = []
user_vectors = []
user_date = []
user_j_vectors = []
for paper in library:
    user_title.append(paper.title)
    vec = paper.vectors.all()[0]
    user_vectors.append(vec.get_vector())
    if paper.journal:
        tmp = paper.journal.vectors.all()[0]
        j_vec = tmp.get_vector()
    else:
        j_vec = np.zeros((128, ))
    user_j_vectors.append(j_vec)

sample_size = 50000

ps_library = Paper.objects.values('pk', 'journal_id')\
    .exclude(pk__in=library).order_by('-pk')[:sample_size]
pvs = PaperVectors.objects.filter(paper__in=ps_library).values('paper__pk', 'paper__title', 'vector')

ps_titles = []
ps_vectors = []
for pv in pvs:
      ps_vectors.append(pv['vector'][:128])
      ps_titles.append(pv['paper__title'])


with open(os.path.join(out_path, 'sample_library.csv'), 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(tuple(['title'] + ['v%s' % i for i in range(128)]))
    for i, title in enumerate(ps_titles):
        row = tuple([title] + [str(val) for val in ps_vectors[i]])
        writer.writerow(row)

with open(os.path.join(out_path, 'user_library.csv'), 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(tuple(['title'] + ['v%s' % i for i in range(128)]))
    for i, title in enumerate(user_title):
        row = tuple([title] + [str(val) for val in user_vectors[i]])
        writer.writerow(row)
