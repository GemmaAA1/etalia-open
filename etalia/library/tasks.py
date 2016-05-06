# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.utils import timezone
from django.conf import settings

from config.celery import celery_app as app
from celery.exceptions import SoftTimeLimitExceeded

from etalia.nlp.models import MostSimilar, PaperNeighbors
from .models import Stats, Paper


@app.task()
def update_stats():
    """Create a new line of stats for the library"""
    stats = Stats.objects.create()
    stats.update()


def get_neighbors_papers(paper_pk, time_span):
    # Get active MostSimilar
    ms = MostSimilar.objects.filter(is_active=True)[0]

    # Get stored neighbors matches
    try:
        neigh_data = PaperNeighbors.objects.get(ms=ms,
                                                paper_id=paper_pk,
                                                time_lapse=time_span)

        if not neigh_data.neighbors or not max(neigh_data.neighbors):
            # neighbors can be filled with zero if none was found previously
            neigh_data.delete()
            raise PaperNeighbors.DoesNotExist
        elif neigh_data.modified > \
                (timezone.now() - timezone.timedelta(
                    days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
            neighbors = neigh_data.neighbors
        else:
            raise PaperNeighbors.DoesNotExist
    except PaperNeighbors.DoesNotExist:   # refresh
        try:
            ms_task = app.tasks['etalia.nlp.tasks.mostsimilar_{name}'.format(
                name=ms.model.name)]
            res = ms_task.apply_async(args=('populate_neighbors',),
                                      kwargs={'paper_pk': paper_pk,
                                              'time_lapse': time_span},
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