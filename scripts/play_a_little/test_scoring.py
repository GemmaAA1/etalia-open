# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from paperstream.users.models import UserLibPaper
from django.contrib.auth import get_user_model
from paperstream.nlp.models import PaperVectors, JournalVectors, Model
import numpy as np
import datetime

model_name = 'dbow-128'
model_name = 'dbow-128-with-words'
id_ex = 43566
id_ex = 363313
id_ex = 404717

model = Model.objects.get(name=model_name)
User = get_user_model()
user = User.objects.first()
feed = user.feeds.first()


seeds = UserLibPaper.objects.filter(userlib=user.lib).order_by('date_created').values('paper_id', 'date_created')
data = PaperVectors.objects.filter(paper__pk__in=seeds.values('paper_id'), model=model).values('paper_id', 'vector')
data_pk = data.values_list('paper_id', flat=True)

data_dic = dict([(d['paper_id'], d['vector'][:model.size]) for d in data])
seeds_l = [(seed['paper_id'], seed['date_created']) for seed in seeds if seed['paper_id'] in data_pk]


seed_mat = np.zeros((len(seeds_l), 128))
date = []
for i, seed in enumerate(seeds_l):
    seed_mat[i, :] = data_dic[seed[0]]
    date.append(seed[1])

# pick target
pv = PaperVectors.objects.get(paper=id_ex, model=model)
vec = np.array(pv.vector[:128])

# compute distance
dist = np.dot(seed_mat, vec)

# save
np.savetxt('dist{id}_{model}.npy'.format(id=id_ex, model=model_name), dist)
with open('date.txt', 'w') as file:
    for item in date:
        file.write("%s\n" % item)

# load

import matplotlib.pyplot as plt
import numpy as np
import datetime

dist = np.loadtxt('dist{id}_{model}.npy'.format(id=id_ex, model=model_name))
date = []
with open('date.txt', 'r') as file:
    while True:
        try:
            date.append(datetime.datetime.strptime(file.readline().strip(), "%Y-%m-%d"))
        except Exception:
            break

import matplotlib.pyplot as plt
plt.plot(date, dist1, '.')
plt.plot(date, dist0, 'r.')
plt.show()

