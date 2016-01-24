# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.utils import timezone
from django.db.models import Q
from django.conf import settings

from paperstream.library.models import Paper
from config.celery import celery_app as app

from .models import AltmetricModel


@app.task()
def update_altmetric_periodic():
    """Update altmetric"""

    # date of paper to update
    d = timezone.datetime.now().date() - timezone.timedelta(days=60)

    # Fetch paper
    ps_pks = Paper.objects\
        .filter(Q(date_ep__gt=d) | (Q(date_ep=None) & Q(date_fs__gt=d)))\
        .values_list('pk', flat=True)[:settings.ALTMETRIC_MAX_PAPERS_PER_PERIOD]

    for pk in ps_pks:
        update_altmetric.apply_async(args=(pk,))


@app.task(rate_limit=1)
def update_altmetric(paper_pk):
    """Celery task for altmetric update"""
    altmetric, _ = AltmetricModel.objects.get_or_create(paper_id=paper_pk)
    altmetric.update()
