# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app
from .models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier, \
    ConsumerJournal, ConsumerPubPeer, ConsumerBiorxiv, ConsumerSpringer, \
    Consumer
from etalia.library.models import Paper
from etalia.core.managers import PaperManager


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
def springer_run_all():
    cps = ConsumerSpringer.objects.all()
    names = [cp.name for cp in cps]
    for name in names:
        springer_run(name)


@app.task()
def biorxiv_run_all():
    cps = ConsumerBiorxiv.objects.all()
    names = [cp.name for cp in cps]
    for name in names:
        biorxiv_run(name)


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

@app.task()
def springer_run(name):
    springer_consumer = ConsumerSpringer.objects.get(name=name)
    springer_consumer.run_once_per_period()

@app.task()
def biorxiv_run(name):
    biorxiv_consumer = ConsumerBiorxiv.objects.get(name=name)
    biorxiv_consumer.run_once_per_period()


@app.task(bind=True, max_retries=3, default_retry_delay=10 * 60)
def populate_journal(self, consumer_id, journal_pk):

    # Find all consumer classes with their type as key
    ConsumerClasses = Consumer.__subclasses__()
    ConsumerClassesType = {}
    for cc in ConsumerClasses:
        ConsumerClassesType[cc.TYPE] = cc

    # Find corresponding Consumer class
    try:
        type_ = Consumer.objects.get(id=consumer_id).type
        consumer = ConsumerClassesType[type_].objects.get(id=consumer_id)
    except Consumer.DoesNotExist:
        raise

    # Populate journal
    try:
        consumer.populate_journal(journal_pk)
    except Exception as exc:
        cj = ConsumerJournal.objects.get(consumer=consumer,
                                         journal_id=journal_pk)
        cj.status = 'retry'
        cj.save(update_fields=['status'])
        raise self.retry(exc=exc, countdown=1)


@app.task()
def populate_pubpeer():
    ppc = ConsumerPubPeer.objects.first()
    ppc.populate()


@app.task(bind=True, rate_limit='1/s')
def consolidate_paper(self, paper_id, **kwargs):
    pm = PaperManager()
    paper = Paper.objects.get(id=paper_id)
    try:
        pm.consolidate_paper(paper, **kwargs)
    except Exception as exc:
        raise self.retry(exc=exc, countdown=5)


@app.task()
def consolidate_library():
    pks = Paper.objects.filter(is_trusted=False).values_list('pk', flat=True)
    for pk in pks:
        consolidate_paper.delay(pk)


@app.task()
def consumer_add(x, y):
    return x + y
