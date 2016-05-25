#!/usr/bin/env python
"""
./update.py
"""
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
import sys
import copy
import yaml
import random
import argparse
from argparse import RawTextHelpFormatter
import names
from subprocess import call
from random import randrange
from datetime import timedelta, datetime, date

NB_THREADS = 50
NB_POSTS = 200
NB_COMMENTS = 200
NB_USERS = 50
NB_RELATIONSHIPS = 100
NB_THREADUSER = 400
INIT_DATA_FILE = os.path.join('data', 'init_data.json')
MODEL_FILE = 'models.yaml'
MODEL_NAME = 'test'


def update():
    update_papers()
    update_users()
    update_threads()
    update_threadusers()
    update_relationships()
    update_mostsimilar()
    update_streams()


def update_papers():
    print('Updating papers...')
    # Update date randomly
    papers = Paper.objects.all()
    end_date = datetime.now().date()
    start_date = end_date - timedelta(days=90)
    for paper in papers:
        paper.date_ep = random_date(start_date, end_date)
        paper.date_fs = random_date(start_date, end_date)
        paper.date_pp = random_date(start_date, end_date)
        paper.save()

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
        altmetric.save()


def update_users():
    print('Updating users...')
    # Fake N users
    if User.objects.all().count() < NB_USERS:
        papers = Paper.objects.all()
        new_user = NB_USERS - User.objects.all().count()
        fake_names = [(names.get_first_name(), names.get_last_name()) for _ in range(new_user)]
        for fn in fake_names:
            user = User.objects.create_user(email='{0}.{1}@fake.com'.format(fn[0], fn[1]),
                                            first_name=fn[0],
                                            last_name=fn[1])
            # add papers
            nb_papers = random.randint(5, 50)
            for i in range(nb_papers):
                ulp, _ = UserLibPaper.objects.get_or_create(
                    userlib=user.lib,
                    paper=papers[random.randrange(0, papers.count())])
                ulp.date_created = random_date(date(2000, 1, 1), date.today())
                ulp.save()
            # add avatar
            ag = AvatarGenerator()
            filename = ag.generate(128, user.email, 'png')
            f = open(filename, 'rb')
            avatar = Avatar(user=user, primary=True)
            avatar.avatar.save(filename, File(f))


def update_threads():
    print('Updating threads...')
    # Fake thread
    if Thread.objects.all().count() < NB_THREADS:
        fixture = AutoFixture(Thread, overwrite_defaults=True)
        fixture.create(NB_THREADS - Thread.objects.all().count())
    if ThreadPost.objects.all().count() < NB_POSTS:
        fixture = AutoFixture(ThreadPost, overwrite_defaults=True)
        fixture.create(NB_POSTS - ThreadPost.objects.all().count())
    if ThreadComment.objects.all().count() < NB_COMMENTS:
        fixture = AutoFixture(ThreadComment, overwrite_defaults=True)
        fixture.create(NB_COMMENTS - ThreadComment.objects.all().count())

    threads = Thread.objects.all()
    posts = ThreadPost.objects.all()
    comments = ThreadComment.objects.all()

    # fix thread logic type
    for thread in threads:
        if thread.type == 1:
            thread.paper = None
            thread.save()
    # fix thread members + Invites
    for post in posts:
        if post.user not in post.thread.members:
            tu, _ = ThreadUser.objects.get_or_create(user=post.user,
                                                  thread=post.thread)
            tu.join()
            if post.thread.privacy == THREAD_PRIVATE:
                # Add Invite accepted
                ThreadUserInvite.objects.get_or_create(
                    thread=post.thread,
                    from_user=post.thread.user,
                    to_user=post.user,
                    status=THREAD_INVITE_ACCEPTED)
    for comment in comments:
        if comment.user not in comment.post.thread.members:
            tu, _ = ThreadUser.objects.get_or_create(user=comment.user,
                                                  thread=comment.post.thread)
            tu.join()
            if comment.post.thread.privacy == THREAD_PRIVATE:
                # Add Invite accepted
                ThreadUserInvite.objects.get_or_create(
                    thread=comment.post.thread,
                    from_user=comment.post.thread.user,
                    to_user=comment.user,
                    status=THREAD_INVITE_ACCEPTED)

    # Embed threads
    pks = Thread.objects.all().exclude(published_at=None).values_list('pk', flat=True)
    model = Model.objects.load(is_active=True)
    model.infer_threads(pks)


def update_threadusers():
    print('Updating ThreadUser states...')
    # ThreadUser
    if ThreadUser.objects.all().count() < NB_THREADUSER:
        fixture = AutoFixture(ThreadUser, overwrite_defaults=True)
        fixture.create(NB_THREADUSER - ThreadUser.objects.all().count())


def update_relationships():
    print('Updating relationships...')
    # Relationships
    if Relationship.objects.all().count() < NB_RELATIONSHIPS:
        fixture = AutoFixture(Relationship, overwrite_defaults=True)
        fixture.create(NB_RELATIONSHIPS - Relationship.objects.all().count())


def update_streams():
    print('Updating streams...')
    us_pk = User.objects.exclude(email__endswith='fake.com').values_list('pk', flat=True)
    # Update users stream
    for user_pk in us_pk:
        reset_stream(user_pk)

    # Update users trend
    for user_pk in us_pk:
        reset_trend(user_pk)


def update_mostsimilar():
    print('Updating MostSimilar...')
    # update mostsimilar
    if not MostSimilar.objects.filter(is_active=True).exists():
        model = Model.objects.get(is_active=True)
        ms = MostSimilar.objects.create(model=model)
        ms.activate()
    mostsimilar_full_update_all()

    # update mostsimilarthread
    if not MostSimilarThread.objects.filter(is_active=True).exists():
        model = Model.objects.get(is_active=True)
        mst = MostSimilarThread.objects.create(model=model)
        mst.activate()
    mostsimilarthread_full_update_all()


def fetch_path():
    # search for config/ is parent tree directory and setup django
    cpath = os.path.abspath(__file__)
    while cpath and not os.path.isdir('config'):
        os.chdir('..')
    sys.path.append(os.path.abspath(os.getcwd()))

    return os.path.abspath(os.getcwd())


def random_date(start, end):
    """
    This function will return a random datetime between two datetime
    objects.
    """
    delta = end - start
    int_delta = (delta.days * 24 * 60 * 60) + delta.seconds
    random_second = randrange(int_delta)
    return start + timedelta(seconds=random_second)


def flush():
    print('Flushing DB...')
    call([os.path.join(root_path, "manage.py"), "flush"])
    sys.exit()


def load():
    print("Loading data from {file}...".format(file=INIT_DATA_FILE))
    call([os.path.join(root_path, "manage.py"), "loaddata",
          os.path.join(root_path, "scripts", "routines", INIT_DATA_FILE)])
    sys.exit()


def fetch_new_papers():
    print('Fetching new papers...')
    pubmed_run('pubmed_all')
    arxiv_run('arxiv_all')
    elsevier_run('elsevier_all')
    sys.exit()


def dump_data():
    print('Dump data to file...')
    call([os.path.join(root_path, "manage.py"), "dumpdata",
          "--exclude=auth",
          "--exclude=contentypes",
          "-o", os.path.join(root_path, "scripts", "routines", INIT_DATA_FILE)])
    sys.exit()


def init_publishers():
    print('Init publishers...')
    call([os.path.join(root_path, "manage.py"), "populate", "publisher", "all"])


def init_consumers():
    print('Init consumers...')
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "pubmed", "--name", "pubmed_all", "--local"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "arxiv", "--name", "arxiv_all"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "elsevier", "--name", "elsevier_all"])


def init_journals():
    print('Init journals...')
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "thomson_local"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "pubmed_local"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "arxiv_local"])


def init_publishers_prod():
    init_publishers()


def init_journals_prod():
    print('Init journals...')
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "thomson"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "pubmed"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "arxiv"])


def init_consumers_prod():
    print('Init consumers...')
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "pubmed", "--name", "pubmed_all"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "arxiv", "--name", "arxiv_all"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "elsevier", "--name", "elsevier_all"])


def init_models():
    print('Init NLP models...')
    # Build NLP model and activate
    with open(os.path.join(root_path, 'scripts', 'routines', MODEL_FILE)) as stream:
        MODELS = yaml.load(stream)
    model_params = [m for m in MODELS if m['name'] == MODEL_NAME][0]

    model = Model.objects.create(**model_params)
    papers = Paper.objects.filter(is_trusted=True,
                                  title__regex=r'^.{5}.*',
                                  abstract__regex=r'^.{10}.*')
    model.dump(papers)
    model.build_vocab_and_train()
    model.activate()

    # Populate paper with NLP signatures
    model.save_journal_vec_from_bulk()
    model.save_paper_vec_from_bulk()


if __name__ == '__main__':

    # Input parser
    parser = argparse.ArgumentParser(description=
                                     'Init and Update Etalia database and related-objects in devolpment:\n',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--init",
                        help="Init database from scratch",
                        action="store_true")
    parser.add_argument("--init-production",
                        help="Init database in production from scratch",
                        action="store_true")
    parser.add_argument("-p", "--papers",
                        help="Fetch new papers",
                        action="store_true")
    parser.add_argument("-f", "--flush",
                        help="Flush database",
                        action="store_true")
    parser.add_argument("-l", "--load",
                        help="Load database with fixture (from {0})".format(INIT_DATA_FILE),
                        action="store_true")
    parser.add_argument("-m", "--models",
                        help="Update NLP models",
                        action="store_true")
    parser.add_argument("-d", "--dump",
                        help="Dump database to file {0}".format(INIT_DATA_FILE),
                        action="store_true")


    args = parser.parse_args()

    # fetch path
    root_path = fetch_path()

    # Reset DB
    if args.flush:
        flush()

    # Load init_data.json
    if args.load:
        load()

    # Dump to init_data.json
    if args.load:
        dump_data()

    # Run pip requirements
    call(["pip",
          "install",
          "-r", os.path.join(root_path, "requirements/development.txt")])

    import django
    django.setup()
    from django.core.files import File
    from django.contrib.auth import get_user_model
    from etalia.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
    from etalia.library.models import Paper
    from etalia.altmetric.tasks import update_altmetric
    from etalia.altmetric.models import AltmetricModel
    from etalia.feeds.tasks import reset_stream, reset_trend
    from etalia.nlp.tasks import mostsimilar_full_update_all
    from etalia.nlp.models import MostSimilarThread, MostSimilar, Model
    from etalia.threads.tasks import mostsimilarthread_full_update_all, \
        embed_threads
    from etalia.threads.models import Thread, ThreadPost, ThreadComment, \
        ThreadUser, ThreadUserInvite
    from etalia.threads.constant import THREAD_PRIVATE, THREAD_INVITE_ACCEPTED, \
        THREAD_PUBLIC
    from etalia.users.models import UserLibPaper, Relationship
    from avatar.models import Avatar
    from utils.avatar import AvatarGenerator
    from autofixture import AutoFixture

    User = get_user_model()

    # Apply migrations
    call([os.path.join(root_path, "manage.py"), "migrate"])

    # Initialization
    if args.init:
        init_publishers()
        init_journals()
        init_consumers()
        fetch_new_papers()
        init_models()

    if args.models:
        init_models()

    if args.init_production:
        init_publishers_prod()
        init_journals_prod()
        init_consumers_prod()

    # Fetch new papers
    if args.papers:
        fetch_new_papers()

    # Fetch data and update etalia objects
    update()
