# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app
from celery.exceptions import SoftTimeLimitExceeded
from time import sleep

from django.utils import timezone
from django.conf import settings
from etalia.nlp.models import PaperEngine, PaperNeighbors
from .models import Stats, Paper


import logging

logger = logging.getLogger(__name__)


@app.task()
def update_stats():
    """Create a new line of stats for the library"""
    stats = Stats.objects.create()
    stats.update()


@app.task(rate_limit='2/s')
def consolidate(paper_pk):
    paper = Paper.objects.get(id=paper_pk)
    paper.consolidate()
    return paper_pk


@app.task()
def consolidate_papers(pks):
    for pk in pks:
        paper = Paper.objects.get(id=pk)
        paper.consolidate_async()


@app.task()
def consolidate_library():
    pks = Paper.objects.filter(is_trusted=False).values_list('pk', flat=True)
    for pk in pks:
        consolidate.delay(pk)


def embed_paper(paper_pk):
    """Send task to embed paper
    """
    from etalia.nlp.tasks import nlp_dispatcher
    nlp_dispatcher.delay('infer_paper', paper_pk)


def embed_papers(pks, model_name, batch_size=1000):
    from etalia.nlp.tasks import nlp_dispatcher
    pks = list(pks)
    nb_papers = len(pks)
    nb_batches = nb_papers // batch_size
    pks_batched = [pks[i*batch_size:(1+i)*batch_size] for i in range(nb_batches)]
    pks_batched.append(pks[nb_batches * batch_size:])

    for i, batch in enumerate(pks_batched):
        nlp_dispatcher.delay('infer_papers', batch)


def get_neighbors_papers(paper_pk, time_span):
    from etalia.nlp.tasks import pe_dispatcher
    # Get active MostSimilar
    pe = PaperEngine.objects.filter(is_active=True)[0]

    # Get stored neighbors matches
    try:
        neigh_data = PaperNeighbors.objects.get(paper_id=paper_pk,
                                                pe=pe,
                                                time_lapse=time_span)
        if not neigh_data.neighbors or not max(neigh_data.neighbors):
            neigh_data.delete()
            raise PaperNeighbors.DoesNotExist
        elif neigh_data.modified > \
                (timezone.now() -
                     timezone.timedelta(days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
            neighbors = neigh_data.neighbors
        else:
            raise PaperNeighbors.DoesNotExist
    except PaperNeighbors.DoesNotExist:   # refresh
        try:
            res = pe_dispatcher.apply_async(args=('populate_neighbors',
                                            paper_pk,
                                            time_span),
                                      timeout=10,
                                      soft_timeout=5)
            neighbors = res.get()
        except SoftTimeLimitExceeded:
            return []
        except KeyError:
            raise

    neigh_pk_list = [neigh for neigh in neighbors[:settings.LIBRARY_NUMBER_OF_NEIGHBORS] if neigh]

    clauses = ' '.join(['WHEN id=%s THEN %s' % (pk, i)
                        for i, pk in enumerate(neigh_pk_list)])
    ordering = 'CASE %s END' % clauses
    return Paper.objects.filter(pk__in=neigh_pk_list).extra(
       select={'ordering': ordering}, order_by=('ordering',))

