# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import yaml
import copy
import random

from paperstream.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
from paperstream.library.models import Paper
from paperstream.altmetric.tasks import update_altmetric
from paperstream.altmetric.models import AltmetricModel
from paperstream.feeds.tasks import reset_stream, reset_trend
from paperstream.nlp.tasks import mostsimilar_update_all

from django.contrib.auth import get_user_model

User = get_user_model()


# POPULATE DATABASE WITH SOME DATA

# fetching some recent papers from the consumers
pubmed_run('pubmed_all')
arxiv_run('arxiv_all')
elsevier_run('elsevier_all')

papers = Paper.objects.filter(is_trusted=True,
                              title__regex=r'^.{5}.*',
                              abstract__regex=r'^.{10}.*')

# Fetch some Altmetric data and fake others (saving time)
alt_objs = []
for paper in papers[:10]:
    update_altmetric(paper.pk)
    obj = copy.copy(paper.altmetric.__dict__)
    obj.pop('paper_id')
    obj.pop('id')
    obj.pop('_state')
    obj.pop('_paper_cache')
    alt_objs.append(obj)
# and randomly assigned to other papers (save time)
for paper in papers:
    alt_obj = alt_objs[random.randint(0, 9)]
    altmetric, _ = AltmetricModel.objects.get_or_create(paper_id=paper.pk)
    for key, value in alt_obj.items():
        setattr(altmetric, key, value)

# update mostsimilar
mostsimilar_update_all()

# Update users stream
us_pk = User.objects.all().values_list('pk', flat=True)
for user_pk in us_pk:
    reset_stream(user_pk)

# Update users trend
for user_pk in us_pk:
    reset_trend(user_pk)

