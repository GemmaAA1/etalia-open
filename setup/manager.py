#!/usr/bin/env python
"""
Script to manage test data and etalia initialization
"""

from __future__ import unicode_literals, absolute_import
import os
import sys
import yaml
import random
import argparse
from argparse import RawTextHelpFormatter
import names
from subprocess import call
from random import randrange
from datetime import timedelta, datetime, date
from django.db import transaction
from django.core.files import File
from django.db.models import Count, F
from autofixture import AutoFixture

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
INIT_DEFAULT_FIXTURE_FILE = os.path.join('data', 'fixture.json')
MODEL_FILE = 'models.yaml'
MODEL_NAME = 'test'


def update():
    """Update etalia (time sensitivity, engines, streams, etc.)"""
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

    # Altmetric score
    for paper in papers:
        altmetric, _ = AltmetricModel.objects.get_or_create(paper_id=paper.pk)
        altmetric.score = float(random.randrange(1, 1000))
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

    # Update oauth user
    users = User.objects.filter(type=USER_INDIVIDUAL).exclude(email__endswith='fake.com')
    for user in users:
        update_oauth_user(user.email)


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


def get_true_users():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .exclude(email__endswith='fake.com')\
        .values_list('id', flat=True)
    return us_pk


def update_streams():
    print('Updating streams...')
    us_pk = get_true_users()
    # Update users stream
    for user_pk in us_pk:
        reset_stream(user_pk)


def update_trends():
    print('Updating trends...')
    us_pk = get_true_users()
    # Update users trend
    for user_pk in us_pk:
        reset_trend(user_pk)


def update_threadfeeds():
    print('Updating threadfeeds...')
    us_pk = get_true_users()
    # Update users threadfeed
    for user_pk in us_pk:
        update_threadfeed(user_pk)


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
          os.path.join(root_path, "setup", file)])


def fetch_new_papers():
    print('Fetching new papers...')
    from etalia.consumers.tasks import populate_journal
    cjs = ConsumerJournal.objects.all()
    for cj in cjs:
        try:
            populate_journal(cj.consumer_id, cj.journal_id)
        except Exception:
            pass


def fetch_new_threads():
    print('Fetching new threads...')
    consumer = ConsumerPubPeer.objects.first()
    consumer.populate()


def dump_data():
    file = os.path.join(root_path, "setup", INIT_DEFAULT_FIXTURE_FILE)
    print('Dump data to {0}'.format(file))
    call([os.path.join(root_path, "manage.py"), "dumpdata",
          "--exclude", "auth",
          "--exclude", "contenttypes",
          "--exclude", "usersession",
          "-o", file])


def init_publishers():
    print('Init publishers...')
    call([os.path.join(root_path, "manage.py"), "populate", "publisher", "all"])


def init_consumers():
    print('Init consumers...')
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "pubmed", "--name", "pubmed_all", "--local"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "arxiv", "--name", "arxiv_all"])
    call([os.path.join(root_path, "manage.py"), "populate", "consumer", "elsevier", "--name", "elsevier_all"])
    ConsumerPubPeer.objects.create()
    c, _ = ConsumerBiorxiv.objects.get_or_create(name='biorxiv')
    j = Journal.objects.get(id_oth='biorxiv')
    cj, _ = ConsumerJournal.objects.get_or_create(consumer=c, journal=j)
    cj.activate()


def init_journals():
    print('Init journals...')
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "thomson_local"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "pubmed_local"])
    call([os.path.join(root_path, "manage.py"), "populate", "journal", "arxiv_local"])
    Journal.objects.get_or_create(title='BioRxiv', id_oth='biorxiv')


def init_models():
    print('Init NLP models...')
    # Build NLP model and activate
    with open(os.path.join(root_path, 'setup', MODEL_FILE)) as stream:
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


def init():
    init_publishers()
    init_journals()
    init_consumers()
    fetch_new_papers()
    init_models()


if __name__ == '__main__':

    # Add path
    root_path = fetch_path()

    import django
    django.setup()
    from django.contrib.auth import get_user_model
    from etalia.library.models import Paper, Journal
    from etalia.altmetric.models import AltmetricModel
    from etalia.feeds.tasks import reset_stream, reset_trend, update_threadfeed
    from etalia.nlp.tasks import paperengine_full_update_all, \
        threadengine_update_all
    from etalia.nlp.models import ThreadEngine, PaperEngine, Model
    from etalia.threads.models import Thread, ThreadPost, ThreadComment, \
        ThreadUser, ThreadUserInvite
    from etalia.threads.constants import THREAD_PRIVATE, \
        THREAD_INVITE_ACCEPTED, THREAD_QUESTION
    from etalia.consumers.models import ConsumerJournal, ConsumerPubPeer, \
        ConsumerBiorxiv
    from etalia.users.models import UserLibPaper, Relationship
    from etalia.users.constants import USER_INDIVIDUAL
    from avatar.models import Avatar
    from setup.utils.avatar import AvatarGenerator

    # Input parser
    parser = argparse.ArgumentParser(
        description=
        'Manage Etalia data and initialization\n'
        '\n'
        'Use cases:\n'
        '  - From a new install WITH fixture, run:\n'
        '      manager.py --load\n'
        '  - To update etalia, run:\n'
        '      manager.py --update\n'
        '  - From a fresh install WITH NO fixture in /data, run:\n'
        '      manager.py --init\n',
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("-i", "--init",
                        help="populate database with test data",
                        action="store_true")
    parser.add_argument("-l", "--load",
                        help="load database from fixture file (default ./{0})"
                        .format(INIT_DEFAULT_FIXTURE_FILE),
                        metavar='file',
                        nargs='?',
                        const=INIT_DEFAULT_FIXTURE_FILE,
                        type=str)
    parser.add_argument("-p", "--papers",
                        help="fetch new papers only",
                        action="store_true")
    parser.add_argument("-u", "--update",
                        help="Update time sensitive data",
                        action="store_true")
    parser.add_argument("-m", "--models",
                        help="init NLP models only",
                        action="store_true")
    parser.add_argument("--user", metavar='email',
                        help="populate OAuth user with data only",
                        type=str)
    parser.add_argument("-d", "--dump",
                        help="dump database to fixture file (default ./{0})"
                        .format(INIT_DEFAULT_FIXTURE_FILE),
                        nargs='?',
                        const=INIT_DEFAULT_FIXTURE_FILE,
                        metavar='file',
                        type=str)
    parser.add_argument("-f", "--flush",
                        help="flush database",
                        action="store_true")
    parser.add_argument("-e", "--engines",
                        help="update engines",
                        action="store_true")
    # Arguments that trigger a standalone action
    UNITARY_UPDATE_ARGS = ['papers',
                           'models',
                           'user',
                           'flush',
                           'dump',
                           'load']

    args = parser.parse_args()
    is_unitary = any([True for arg in UNITARY_UPDATE_ARGS
                      if not getattr(args, arg) in [None, False]])

    User = get_user_model()

    # SINGLE UPDATE
    if args.user:
        update_oauth_user(args.user)

    if args.models:
        init_models()

    # Reset DB
    if args.flush:
        flush()

    # Load init_data.json
    if args.load:
        load(args.load)
        update()

    # Dump to init_data.json
    if args.dump:
        dump_data()

    # Fetch new papers
    if args.papers:
        fetch_new_papers()

    # Update engines
    if args.engines:
        update_engines()

    # Initialization
    if args.init:
        init()
        update()

    # Update
    if args.update:
        if not Paper.objects.count():
            default_fixture = os.path.join(root_path, "setup",
                                           INIT_DEFAULT_FIXTURE_FILE)
            if os.path.isfile(default_fixture):
                load(INIT_DEFAULT_FIXTURE_FILE)
            else:
                init()
        update()
