# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from paperstream.users.models import UserLibPaper
from django.contrib.auth import get_user_model
from paperstream.nlp.models import PaperVectors, JournalVectors, Model
from paperstream.library.models import Paper
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.cluster import AffinityPropagation, DBSCAN
import numpy as np
import os
import datetime

model_name = 'dbow-128'
model_name = 'dbow-128-with-words'
model_name = 'dm-128'

save_path = '/home/ubuntu/staging/test'
id_exs = [862220, 413285, 43566, 566276, 531410]
titles = ['Dietary and lifestyle advice for pregnant women who are overweight or obese: the LIMIT randomized trial.',
 'Cluster analysis of quantitative parametric maps from DCE-MRI: application in evaluating heterogeneity of tumor response to antiangiogenic treatment.',
 'Supervised Fine Tuning for Word Embedding with Integrated Knowledge',
 "Tomoregulin (TMEFF2) Binds Alzheimer's Disease Amyloid-beta (Abeta) Oligomer and AbetaPP and Protects Neurons from Abeta-Induced Toxicity.",
 'Prefrontal Engagement and Reduced Default Network Suppression Co-occur and Are Dynamically Coupled in Older Adults: The Default-Executive Coupling Hypothesis of Aging.']


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

# for id_ex in id_exs:
# pick target
id_ex = id_exs[4]
pv = PaperVectors.objects.get(paper=id_ex, model=model)
vec = np.array(pv.vector[:128])

# compute distance
dist = np.dot(seed_mat, vec)

# save
np.savetxt(os.path.join(save_path, 'dist{id}_{model}.npy'.format(id=id_ex, model=model_name)), dist)
with open('date.txt', 'w') as file:
    for item in date:
        file.write("%s\n" % item)


# get article title
ps = dict(Paper.objects.filter(pk__in=id_exs).values_list('pk', 'title'))
titles = [ps[pk] for pk in id_exs]


# load

import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import datetime
import seaborn

dists = []
for id_ex in id_exs:
    dists.append(np.loadtxt('dist{id}_{model}.npy'.format(id=id_ex, model=model_name)))
    date = []
    with open('date.txt', 'r') as file:
        while True:
            try:
                date.append(datetime.datetime.strptime(file.readline().strip(), "%Y-%m-%d"))
            except Exception:
                break

plt.ion()
for i, dist in enumerate(dists):
    plt.figure()
    plt.plot(date, dist, '.')
    plt.ylim([-0.7, 0.7])
    plt.title(titles[i])
    plt.show()


day_lapse = np.array([(d.date() - datetime.datetime.now().date()).days for d in date])
baseline = 0.5
k = 0.1
delay = -60.0
weigh = baseline + (1-baseline) / (1 + np.exp(- (day_lapse - delay) * k))

# distsw = []
# for dist in dists:
#     distsw.append(weigh * np.sign(dist) * dist ** 2 * 10)
#
# dists = distsw

# compute histogram
bins = np.arange(-.7, .7, .01)
hist = []
skews = []
stds = []
means = []

for dist in dists:
    means.append(np.mean(dist))
    skews.append(stats.skew(dist))
    stds.append(np.std(dist))
    hist.append(np.histogram(dist, bins))

plt.figure()
col = ['g', 'r', 'b', 'm', 'y']
for i, dist in enumerate(dists):
    plt.figure()
    n, bins, patches = plt.hist(dist, bins, normed=1, histtype='stepfilled')
    plt.setp(patches, 'facecolor', col[i], 'alpha', 0.6)
    plt.title(titles[i])
    plt.show()
plt.legend(['unrelated', 'DCE related', 'word2vec related'])
plt.show()

i = 0
plt.plot(hist[i][1], hist[i][0])
plt.show()


plt.figure()
for i, dist in enumerate(dists[:3]):
    n, bins, patches = plt.hist(dist, bins, normed=1, histtype='stepfilled')
    # plt.setp(patches, 'facecolor', 'alpha', 0.6)
plt.legend(['unrelated', 'DCE related', 'word2vec related'])
plt.show()



from django.contrib.auth import get_user_model
from paperstream.feeds.scoring import ContentBasedSimple
from sklearn import manifold, datasets

import  matplotlib.pylab as plt

User = get_user_model()
us = User.objects.all()
user = us[0]

stream = user.streams.first()
self = ContentBasedSimple(stream=stream)

self.build_profile_ind2jourpk()
self.build_profile_ind2authpk()
date_vec = self.build_date_vec()
# build seed mat
seed_vec_mat = self.build_paper_vec_mat(self.seed_data)

# build seed author mat
seed_auth_mat = self.build_auth_utility_mat(self.seed_auth_data)

# build journal mat
seed_jour_mat = self.build_jour_utility_mat(self.seed_data)

# concatenate these 3 mats
seed_mat = np.hstack((self.vec_w * seed_vec_mat,
                      self.auth_w * seed_auth_mat,
                      self.jour_w * seed_jour_mat))

for i in range(seed_mat.shape[1]):
    tmp = seed_mat[:,i]
    seed_mat[:,i] = (tmp - np.mean(tmp))/np.std(tmp)


X = seed_mat[:, :128]

# normalize
norm = np.linalg.norm(X)
if norm > 0:
    X /= norm


tsne = manifold.TSNE(n_components=2, init='pca', random_state=0,
                     perplexity=50,
                     early_exaggeration=8.0,
                     learning_rate=200,
                     verbose=2,
                     method='exact')
Y = tsne.fit_transform(X)

fig = plt.figure()
color = []
plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
plt.axis('tight')
plt.show()


cl = MiniBatchKMeans(n_clusters=5,
                     init='k-means++',
                     n_init=1,
                     init_size=100,
                     batch_size=100)

cl.fit(Y)

# cl = DBSCAN(eps=0.3*6, min_samples=20).fit(Y)
# cl = AffinityPropagation(preference=-50).fit(Y)

color = cl.labels_
fig = plt.figure()
plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
plt.axis('tight')
plt.show()

