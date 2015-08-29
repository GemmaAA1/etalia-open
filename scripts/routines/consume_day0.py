# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from paperstream.consumers.models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier
from django.db.models import Q
import time


def init_consumer(consumers):

    for consumer in consumers:
        consumer.activate_all()
        # Get journal active
        consumer_journal_active = consumer.consumerjournal_set.filter(
            Q(status='idle') | Q(status='consuming') | Q(status='in_queue'))
        # queue journal for consumption
        err = []
        for consumerjournal in consumer_journal_active:
            try:
                consumer.populate_journal(consumerjournal.journal.pk)
                time.sleep(1)
            except Exception as msg:
                err.append({consumerjournal.journal.pk: msg})

    return err

if __name__ == '__main__':

    pubmed = ConsumerPubmed.objects.first()
    elsevier = ConsumerElsevier.objects.first()
    arxiv = ConsumerArxiv.objects.first()

    cons = [pubmed, elsevier, arxiv]

    init_consumer(cons)