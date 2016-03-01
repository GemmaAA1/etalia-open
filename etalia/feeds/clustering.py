# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import numpy as np
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.cluster import MeanShift, estimate_bandwidth
from django.conf import settings
from sklearn import manifold

from etalia.nlp.models import PaperVectors, Model


class Clustering(object):
    """ Clustering """

    def __init__(self, stream, **kwargs):
        self.stream = stream

        self.is_ready = False
        self.model = Model.objects.get(is_active=True)

        # Init seed data
        self.seed_pks = self.get_seed_pks()
        self.seed_data = self.get_data(self.seed_pks)

    def get_seed_pks(self):
        """Retrieve seed paper pk"""
        #TODO: optimze SQL query

        seed_pks = PaperVectors.objects\
            .filter(paper__pk__in=self.stream.seeds.values('pk'),
                    model=self.stream.user.settings.stream_model)\
            .values_list('paper__pk', flat=True)

        return list(seed_pks)

    def get_data(self, pks):
        """Get paper metadata and vector and sort by pks"""
        data = list(PaperVectors.objects\
            .filter(
                paper__pk__in=pks,
                model=self.model)\
            .values('vector',
                    'paper__pk',
                    'paper__journal__pk'))
        # sort
        data.sort(key=lambda d: pks.index(d['paper__pk']))

        return data

    def build_paper_vec_mat(self, data):
        """Build 2D array of paper vectors (matches x vector)"""
        # init seed-mat
        vectors = [d['vector'][:self.model.size] for d in data]
        return np.array(vectors)

    def cluster(self):

        # Build vec matrix
        seed_vec_mat = self.build_paper_vec_mat(self.seed_data)

        # dimensionality reduction
        tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
        Y = tsne.fit_transform(seed_vec_mat)

        # # Perform clustering
        km = MiniBatchKMeans(n_clusters=settings.FEED_NB_CLUSTERS,
                             init='k-means++',
                             n_init=1,
                             init_size=1000,
                             batch_size=1000)

        km.fit(seed_vec_mat)



        # clustering
        bandwidth = estimate_bandwidth(Y, quantile=0.2, n_samples=500)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(Y)

        return self.seed_pks, ms

