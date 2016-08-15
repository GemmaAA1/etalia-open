#!/usr/bin/env python
"""
./update.py --help for help
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
from django.db import transaction

NB_THREADS = 50
NB_POSTS = 200
NB_COMMENTS = 200
NB_USERS = 50
NB_OAUTH_INVITE = 5
NB_RELATIONSHIPS = 100
NB_THREADUSER = 400
NB_THREADS_USER = 10
NB_POSTS_USER = 20
NB_COMMENTS_USER = 20
NB_RELATIONSHIPS_USER = 20
NB_THREADUSER_USER = 20
INIT_DEFAULT_FIXTURE_FILE = os.path.join('data', 'init_data.json')
MODEL_FILE = 'models.yaml'
MODEL_NAME = 'test'


def update():
    update_papers()
    update_users()
    update_threads()
    update_threadusers()
    update_relationships()
    update_engines()
    update_streams()
    update_trends()
    update_threadfeeds()


def update_oauth_user(email):
    print('Updating user ({email})'.format(email=email))
    user = User.objects.get(email=email)

    # Thread
    fixture = AutoFixture(Thread,
                          field_values={'user': user},
                          overwrite_defaults=True)
    if Thread.objects.filter(user=user).count() < NB_THREADS_USER:
        fixture.create(NB_THREADS_USER -
                       Thread.objects.filter(user=user).count())
    # Posts
    fixture = AutoFixture(ThreadPost,
                          field_values={'user': user},
                          overwrite_defaults=True)
    if ThreadPost.objects.filter(user=user).count() < NB_POSTS_USER:
        fixture.create(NB_POSTS_USER -
                       ThreadPost.objects.filter(user=user).count())
    # Comments
    fixture = AutoFixture(ThreadComment,
                          field_values={'user': user},
                          overwrite_defaults=True)
    if ThreadComment.objects.filter(user=user).count() < NB_COMMENTS_USER:
        fixture.create(NB_COMMENTS_USER -
                       ThreadComment.objects.filter(user=user).count())

    threads = Thread.objects.filter(user=user)
    posts = ThreadPost.objects.filter(user=user)
    comments = ThreadComment.objects.filter(user=user)

    # Embed threads
    pks = threads.exclude(published_at=None).values_list('pk', flat=True)
    model = Model.objects.load(is_active=True)
    model.infer_threads(pks)

    fix_thread_logic(threads, posts, comments)

    # ThreadUser
    if ThreadUser.objects.filter(user=user).count() < NB_THREADUSER_USER:
        fixture = AutoFixture(ThreadUser, overwrite_defaults=True,
                              field_values={'user': user})
        fixture.create(NB_THREADUSER_USER -
                       ThreadUser.objects.filter(user=user).count())

    # Add invites
    if ThreadUserInvite.objects.filter(to_user=user).count() < NB_OAUTH_INVITE:
        nb_threads = Thread.objects.all().exclude(user=user).count()
        idx = [random.randint(0, nb_threads-1) for _ in range(NB_OAUTH_INVITE)]
        threads = [Thread.objects.all().exclude(user=user)[i] for i in idx]
        for thread in threads:
            ThreadUserInvite.objects.get_or_create(
                        thread=thread,
                        from_user=thread.user,
                        to_user=user)

    # Relationships
    if Relationship.objects.filter(from_user=user).count() < NB_RELATIONSHIPS_USER:
        fixture = AutoFixture(Relationship, overwrite_defaults=True)
        fixture.create(NB_RELATIONSHIPS -
                       Relationship.objects.filter(from_user=user).count())

    # Update streams
    reset_stream(user.pk)
    reset_trend(user.pk)

    sys.exit()


def update_papers():
    print('Updating papers...')
    # Update date randomly
    with transaction.atomic():
        papers = Paper.objects.select_for_update().all()
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

    # Embed threads
    pks = threads.exclude(published_at=None).values_list('pk', flat=True)
    model = Model.objects.load(is_active=True)
    model.infer_threads(pks)

    fix_thread_logic(threads, posts, comments)


def fix_thread_logic(threads, posts, comments):

    # fix thread logic type
    for thread in threads:
        if thread.type == THREAD_QUESTION:
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


def update_trends():
    print('Updating trends...')
    us_pk = User.objects.exclude(email__endswith='fake.com').values_list('pk', flat=True)
    # Update users trend
    for user_pk in us_pk:
        reset_trend(user_pk)


def update_threadfeeds():
    print('Updating threadfeeds...')
    us_pk = User.objects.exclude(email__endswith='fake.com').values_list('pk', flat=True)
    # Update users threadfeed
    for user_pk in us_pk:
        reset_threadfeed(user_pk)


def update_engines():
    print('Updating Engines...')
    # update PaperEngine
    if not PaperEngine.objects.filter(is_active=True).exists():
        model = Model.objects.get(is_active=True)
        pe = PaperEngine.objects.create(model=model)
        pe.activate()
    paperengine_full_update_all()

    # update ThreadEngine
    if not ThreadEngine.objects.filter(is_active=True).exists():
        model = Model.objects.get(is_active=True)
        te = ThreadEngine.objects.create(model=model)
        te.activate()
    threadengine_update_all()


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


def load(file):
    print("Loading data from {file}...".format(file=file))
    call([os.path.join(root_path, "manage.py"), "loaddata",
          os.path.join(root_path, "scripts", "routines", file)])


def fetch_new_papers():
    print('Fetching new papers...')
    pubmed_run('pubmed_all')
    arxiv_run('arxiv_all')
    elsevier_run('elsevier_all')


def dump_data():
    print('Dump data to file...')
    call([os.path.join(root_path, "manage.py"), "dumpdata",
          "--exclude", "auth", "--exclude", "contenttypes",
          "-o", os.path.join(root_path, "scripts", "routines", INIT_DEFAULT_FIXTURE_FILE)])


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

    if not Model.objects.filter(name=model_params['name']).exists():

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
    else:
        print('NLP Model already setup!')

if __name__ == '__main__':


    # Input parser
    parser = argparse.ArgumentParser(
        description=
        'Update Etalia DB and Engines in development\n'
        '\n'
        'Use cases:\n'
        '  - When starting from scratch, run:\n'
        '      update.py -l (load test data and update engines)\n'
        '  - When not starting from scratch, run:\n'
        '      update.py (update engines and DB)\n'
        '  - If previous fails, run:\n'
        '      update.py -f (flush the DB)\n'
        '      update.py -l (reload test data and update engines)\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--init",
                        help="Populate database with test data",
                        action="store_true")
    parser.add_argument("-l", "--load",
                        help="Load database from fixture file (default ./{0})".format(INIT_DEFAULT_FIXTURE_FILE),
                        metavar='file',
                        nargs='?',
                        const=INIT_DEFAULT_FIXTURE_FILE,
                        type=str)
    parser.add_argument("-p", "--papers",
                        help="Fetch new papers only",
                        action="store_true")
    parser.add_argument("-u", "--users",
                        help="Populate database with test users",
                        action="store_true")
    parser.add_argument("-m", "--models",
                        help="Init NLP models only",
                        action="store_true")
    parser.add_argument("--user", metavar='email',
                        help="Populate OAuth user with data only",
                        type=str)
    parser.add_argument("-d", "--dump",
                        help="Dump database to fixture file (default ./{0})".format(INIT_DEFAULT_FIXTURE_FILE),
                        nargs='?',
                        const=INIT_DEFAULT_FIXTURE_FILE,
                        metavar='file',
                        type=str)
    parser.add_argument("-f", "--flush",
                        help="Flush database",
                        action="store_true")
    parser.add_argument("-e", "--engines",
                        help="Update engines",
                        action="store_true")
    parser.add_argument("--init-production",
                        help="Init database in production from scratch",
                        action="store_true")
    UNITARY_UPDATE_ARGS = ['papers', 'models', 'user', 'flush', 'dump', 'load']

    args = parser.parse_args()
    is_unitary = any([True for arg in UNITARY_UPDATE_ARGS
                      if not getattr(args, arg) in [None, False]])

    # fetch path
    root_path = fetch_path()

    # Run pip requirements
    if not is_unitary:
        call(["pip",
              "install",
              "-r", os.path.join(root_path, "requirements/development.txt")])

    # Apply migrations
    call([os.path.join(root_path, "manage.py"), "migrate"])

    import django
    django.setup()
    from django.core.files import File
    from django.contrib.auth import get_user_model
    from etalia.consumers.tasks import pubmed_run, arxiv_run, elsevier_run
    from etalia.library.models import Paper
    from etalia.altmetric.tasks import update_altmetric
    from etalia.altmetric.models import AltmetricModel
    from etalia.feeds.tasks import reset_stream, reset_trend, reset_threadfeed
    from etalia.nlp.tasks import paperengine_full_update_all, threadengine_update_all
    from etalia.nlp.models import ThreadEngine, PaperEngine, Model
    from etalia.threads.models import Thread, ThreadPost, ThreadComment, \
        ThreadUser, ThreadUserInvite
    from etalia.threads.constant import THREAD_PRIVATE, THREAD_INVITE_ACCEPTED, \
        THREAD_PAPER, THREAD_QUESTION
    from etalia.users.models import UserLibPaper, Relationship
    from avatar.models import Avatar
    from utils.avatar import AvatarGenerator
    from autofixture import AutoFixture

    User = get_user_model()

    # SINGLE UPDATE
    if args.user:
        update_oauth_user(args.user)
        sys.exit()

    if args.models:
        init_models()
        sys.exit()

    # Reset DB
    if args.flush:
        flush()
        sys.exit()

    # Load init_data.json
    if args.load:
        load(args.load)
        sys.exit()

    # Dump to init_data.json
    if args.dump:
        dump_data()
        sys.exit()

    # Fetch new papers
    if args.papers:
        fetch_new_papers()
        sys.exit()

    # Update engines
    if args.engines:
        update_engines()
        sys.exit()

    # MASTER UPDATE
    # Initialization in local environment
    if args.init:
        init_publishers()
        init_journals()
        init_consumers()
        fetch_new_papers()
        init_models()

    if args.init_production:
        init_publishers_prod()
        init_journals_prod()
        init_consumers_prod()

    # Fetch data and update etalia objects
    update()
