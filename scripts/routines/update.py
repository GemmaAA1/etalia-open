#!/usr/bin/env python
"""
./update.py
"""
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
import sys
import copy
import random
from subprocess import call


def update():
    from django.contrib.auth import get_user_model
    from etalia.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
    from etalia.library.models import Paper
    from etalia.altmetric.tasks import update_altmetric
    from etalia.altmetric.models import AltmetricModel
    from etalia.feeds.tasks import reset_stream, reset_trend
    from etalia.nlp.tasks import mostsimilar_full_update_all

    User = get_user_model()

    # POPULATE DATABASE WITH SOME DATA

    # fetching some recent papers from the consumers
    pubmed_run('pubmed_all')
    arxiv_run('arxiv_all')
    elsevier_run('elsevier_all')

    papers = Paper.objects.all()

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
            if key == 'score':
                value = float(random.randrange(1, 1000))
            setattr(altmetric, key, value)

    # update mostsimilar
    mostsimilar_full_update_all()

    # Update users stream
    us_pk = User.objects.all().values_list('pk', flat=True)
    for user_pk in us_pk:
        reset_stream(user_pk)

    # Update users trend
    for user_pk in us_pk:
        reset_trend(user_pk)


def fetch_path():
    # search for config/ is parent tree directory and setup django
    cpath = os.path.abspath(__file__)
    while cpath and not os.path.isdir('config'):
        os.chdir('..')
    sys.path.append(os.path.abspath(os.getcwd()))

    return os.path.abspath(os.getcwd())


def setup_django():
    import django
    django.setup()


if __name__ == '__main__':
    """Update db and relevant objects with fresh data
    """

    # fetch path
    root_path = fetch_path()

    # Run pip requirements
    call(["pip",
          "install",
          "-r", os.path.join(root_path, "requirements/development.txt")])

    # Apply migrations
    call([os.path.join(root_path, "manage.py"), "migrate"])

    # Setup django
    setup_django()

    # Fetch data and update etalia objects
    update()
