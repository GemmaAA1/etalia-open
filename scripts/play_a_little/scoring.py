# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import csv
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

data_path = '/Users/nicolaspannetier/Projects/paperstream/data'
user_file = 'nicolas.pannetier@gmail.com_library.csv'
lib_file = 'sample_library.csv'


def normalize(a):
    return a / np.linalg.norm(a)

# Read files
lib_mat = []
lib_titles = []
with open(os.path.join(data_path, lib_file), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    # trash first line
    reader.__next__()
    for row in reader:
        lib_titles.append(row[0])
        lib_mat.append(row[1:])

user_mat = []
user_titles = []
with open(os.path.join(data_path, user_file), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    # trash first line
    reader.__next__()
    for row in reader:
        user_titles.append(row[0])
        user_mat.append(row[1:])


# remove duplicate
lib_titles, ind = np.unique(np.array(lib_titles), return_index=True)
lib_mat = np.array(lib_mat, dtype=np.float)[ind, :]
user_titles, ind = np.unique(np.array(user_titles), return_index=True)
user_mat = np.array(user_mat, dtype=np.float)[ind, :]

# re-normalize
lib_mat = np.apply_along_axis(normalize, 1, lib_mat)
user_mat = np.apply_along_axis(normalize, 1, user_mat)


# Get baseline stats
nb_samples = 10000
samples = np.unique(np.random.randint(0, len(lib_titles), nb_samples))

# Get baseline dist distribution
lib_mat_red = lib_mat[samples, :]
dist_lib = 1.0 - np.dot(lib_mat_red, lib_mat_red.T)
tri = np.tri(dist_lib.shape[0], k=-1, dtype=np.bool)
dist_lib_tri = dist_lib[tri]

step = 0.01
bins = np.arange(0., 2., step)
dist_lib_hist = np.histogram(dist_lib_tri, bins, density=True)

# compare user paper to baseline
dist_user = 1.0 - np.dot(lib_mat_red, user_mat.T)
kl = []
hists = []
for i in range(user_mat.shape[0]):
    d = np.histogram(dist_user[:, i], bins, density=True)
    kl.append(stats.entropy(d[0], dist_lib_hist[0]))
    hists.append(d)


plt.figure()
plt.imshow(user_mat)
plt.show()



# SCORING

dot = np.dot(lib_mat, user_mat.T)

samples = np.random.randint(0, len(lib_titles), 10000)

i = 0
p = dot[samples, i]


# scoring sort
dist = 1.0 - dot
# get n closest neighbors per paper in user lib
cutoff = 5
ind = np.argpartition(dist, cutoff, axis=0)[:cutoff, :][:]
ind_unique = np.unique(ind[:])

# count occurrence
occ = []
for i, idx in enumerate(ind_unique):
    occ.append((idx, np.sum(ind == idx)))
# normalize by number of paper in user lib
occ = list(map(lambda x: (x[0], x[1]/(user_mat.shape[0]*cutoff)), occ))
occ_sorted = sorted(occ, key=lambda x: x[1], reverse=True)
for i, x in enumerate(occ_sorted[:50]):
    print('[{0}]: {2} {1}'.format(x[1], lib_titles[x[0]], i))

# count sum up distance
occ = []
for i, idx in enumerate(ind_unique):
    occ.append((idx, np.mean(dist[idx, np.where(ind == idx)[1]])))
occ_sorted = sorted(occ, key=lambda x: x[1], reverse=False)
for i, x in enumerate(occ_sorted[:50]):
    print('[{0}]: {2} {1}'.format(x[1], lib_titles[x[0]], i))


occ_a = np.array([o[1] for o in occ])



dist = 1.0 - dot
ind = np.where(dist < 0.4)

eigs_v = []
sdp = []

for i in range(100):
    vec1 = lib_mat_red[i, :]
    c = np.outer(vec1, vec1)
    eigs = np.linalg.eigvals(c)
    eigs_v.append(eigs)
    sdp.append(np.all(eigs > 0 ))


#
dist16 = np.dot(user_mat[16,:], lib_mat_red.T)
dist2 = np.dot(user_mat[2,:], lib_mat_red.T)
dist3 = np.dot(user_mat[3,:], lib_mat_red.T)

step = 0.01
bins = np.arange(-1., 1., step)
d16 = np.histogram(dist16, bins, density=True)
d2 = np.histogram(dist2, bins, density=True)
d3 = np.histogram(dist3, bins, density=True)

stats.entropy(d16[0][np.where(d16[0] > 0 and d2[0] > 0)], d2[0])
stats.entropy(d16[0], d3[0])
lab = ['{0:05d}'.format(b) for b in range(len(bins)-1)]
foo = dit.Distribution(lab, d16[0].tolist())
bar = dit.Distribution(['A','B','C'],[0.1,0.0,0.9])
dit.divergences.jensen_shannon_divergence([foo,bar])