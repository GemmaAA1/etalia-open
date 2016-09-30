# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.utils import timezone
from django.conf import settings

from celery.exceptions import SoftTimeLimitExceeded

from etalia.nlp.models import Model, ThreadEngine
from config.celery import celery_app as app
from etalia.nlp.models import ThreadNeighbors
from .models import Thread
from .utils import send_invite_thread_email


import logging

logger = logging.getLogger(__name__)


def embed_thread(thread_pk):
    """Send task to embed thread
    """
    try:
        from etalia.nlp.tasks import nlp_dispatcher
        nlp_dispatcher.delay('infer_thread', thread_pk)
    except ImportError:
        pass


def embed_threads(pks, batch_size=1000):
    from etalia.nlp.tasks import nlp_dispatcher
    pks = list(pks)
    nb_papers = len(pks)
    nb_batches = nb_papers // batch_size
    pks_batched = [pks[i * batch_size:(1 + i) * batch_size] for i in
                   range(nb_batches)]
    pks_batched.append(pks[nb_batches * batch_size:])

    for batch in pks_batched:
        nlp_dispatcher.delay('infer_threads', batch)


def get_neighbors_threads(thread_pk, time_span):
    from etalia.nlp.tasks import te_dispatcher
    # Get active ThreadEngine
    te = ThreadEngine.objects.filter(is_active=True)[0]

    # Get stored neighbors matches
    try:
        neigh_data = ThreadNeighbors.objects.get(thread_id=thread_pk,
                                                 te=te,
                                                 time_lapse=time_span)
        if not neigh_data.neighbors or not max(
                neigh_data.neighbors):  # neighbors can be filled with zero if none was found previously
            neigh_data.delete()
            raise ThreadNeighbors.DoesNotExist
        elif neigh_data.modified > (timezone.now() - timezone.timedelta(
                days=settings.NLP_NEIGHBORS_REFRESH_TIME_LAPSE)):
            neighbors = neigh_data.neighbors
        else:
            raise ThreadNeighbors.DoesNotExist
    except ThreadNeighbors.DoesNotExist:  # refresh
        try:
            res = te_dispatcher.apply_async(args=('populate_neighbors',
                                            thread_pk,
                                            time_span),
                                      timeout=10,
                                      soft_timeout=5)
            neighbors = res.get()
        except SoftTimeLimitExceeded:
            neighbors = []
        except KeyError:
            raise

    neigh_pk_list = [neigh for neigh in neighbors[:settings.THREADS_NUMBER_OF_NEIGHBORS]
                     if neigh]

    clauses = ' '.join(['WHEN id=%s THEN %s' % (pk, i)
                        for i, pk in enumerate(neigh_pk_list)])
    ordering = 'CASE %s END' % clauses
    return Thread.objects.filter(pk__in=neigh_pk_list).extra(
        select={'ordering': ordering}, order_by=('ordering',))


@app.task()
def async_send_invite_thread_email(from_id, to_id, thread_id):
    """Send email for thread invite"""
    send_invite_thread_email(from_id, to_id, thread_id)
