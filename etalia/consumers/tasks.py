# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app
from .models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier, \
    ConsumerJournal


@app.task()
def pubmed_run_all():
    cps = ConsumerPubmed.objects.all()
    names = [cp.name for cp in cps]
    for name in names:
        pubmed_run(name)


@app.task()
def arxiv_run_all():
    cps = ConsumerArxiv.objects.all()
    names = [cp.name for cp in cps]
    for name in names:
        arxiv_run(name)


@app.task()
def elsevier_run_all():
    cps = ConsumerElsevier.objects.all()
    names = [cp.name for cp in cps]
    for name in names:
        elsevier_run(name)


@app.task()
def pubmed_run(name):
    pubmed_consumer = ConsumerPubmed.objects.get(name=name)
    pubmed_consumer.run_once_per_period()


@app.task()
def arxiv_run(name):
    arxiv_consumer = ConsumerArxiv.objects.get(name=name)
    arxiv_consumer.run_once_per_period()


@app.task()
def elsevier_run(name):
    elsevier_consumer = ConsumerElsevier.objects.get(name=name)
    elsevier_consumer.run_once_per_period()


@app.task(bind=True)
def populate_journal(self, consumer_id, type, journal_pk):
    try:
        if type == 'PUB':
            ConsumerClass = ConsumerPubmed
        elif type == 'ELS':
            ConsumerClass = ConsumerElsevier
        elif type == 'ARX':
            ConsumerClass = ConsumerArxiv
        else:
            raise ValueError('Consumer type unknown {0}'.format(type))
        consumer = ConsumerClass.objects.get(id=consumer_id)
        consumer.populate_journal(journal_pk)
    except Exception as exc:
        cj = ConsumerJournal.objects.get(consumer=consumer,
                                         journal_id=journal_pk)
        cj.status = 'retry'
        cj.save(update_fields=['status'])
        raise self.retry(exc=exc, countdown=1)


@app.task()
def consumer_add(x, y):
    return x + y
