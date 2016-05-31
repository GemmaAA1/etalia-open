# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app

from celery.exceptions import SoftTimeLimitExceeded

from django.utils import timezone
from django.conf import settings
from .models import Stats, Paper
from etalia.nlp.models import Model, PaperEngine, PaperNeighbors

import logging

logger = logging.getLogger(__name__)


@app.task()
def update_stats():
    """Create a new line of stats for the library"""
    stats = Stats.objects.create()
    stats.update()


def embed_paper(paper_pk):
    """Send chain task to embed paper
    """
    model_names = Model.objects\
        .filter(is_active=True)\
        .values_list('name', flat=True)
    for model_name in model_names:
        # Send task for embedding
        try:
            model_task = app.tasks['etalia.nlp.tasks.pe_{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Model task for {model_name} not defined'.format(
                model_name=model_name))
            continue

        model_task.delay('infer_paper', paper_pk)


def embed_papers(pks, model_name, batch_size=1000):
    try:
        model_task = app.tasks['etalia.nlp.tasks.pe_{model_name}'.format(
            model_name=model_name)]
    except KeyError:
        logger.error('Embeding task for {model_name} not defined'.format(
            model_name=model_name))
        raise KeyError
    pks = list(pks)
    nb_papers = len(pks)
    nb_batches = nb_papers // batch_size
    pks_batched = [pks[i*batch_size:(1+i)*batch_size] for i in range(nb_batches)]
    pks_batched.append(pks[nb_batches * batch_size:])

    for batch in pks_batched:
        model_task.delay('infer_papers', batch)


def get_neighbors_papers(paper_pk, time_span):

    # Get active MostSimilar
    pe = PaperEngine.objects.filter(is_active=True)

    # Get stored neighbors matches
    try:
        neigh_data = PaperNeighbors.object.get(paper_id=paper_pk,
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
            pe_task = app.tasks[
                'etalia.nlp.tasks.pe_{name}'.format(name=pe.model.name)]
            res = pe_task.apply_async(args=('populate_neighbors',
                                            paper_pk,
                                            time_span),
                                      timeout=10,
                                      soft_timeout=5)
            neighbors = res.get()
        except SoftTimeLimitExceeded:
            neighbors = []
        except KeyError:
            raise

    neigh_pk_list = [neigh for neigh in neighbors[:settings.NUMBER_OF_NEIGHBORS] if neigh]

    clauses = ' '.join(['WHEN id=%s THEN %s' % (pk, i)
                        for i, pk in enumerate(neigh_pk_list)])
    ordering = 'CASE %s END' % clauses
    return Paper.objects.filter(pk__in=neigh_pk_list).extra(
       select={'ordering': ordering}, order_by=('ordering',))
