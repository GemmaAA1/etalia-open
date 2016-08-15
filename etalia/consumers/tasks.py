# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from config.celery import celery_app as app
from .models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier
from etalia.library.models import Journal
from django.conf import settings


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
def populate_journal(self, consumer, journal_pk):
    try:
        return consumer.populate_journal(journal_pk)
    except Exception as exc:
        raise self.retry(exc=exc, countdown=60)
