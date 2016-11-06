# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from django.db import transaction
from config.celery import celery_app as app
from django.db.models import F, Count
from django.contrib.auth import get_user_model
from django.utils import timezone
from etalia.library.models import Paper
from etalia.threads.models import Thread
from etalia.last_seen.models import LastSeen
from .models import Stream, Trend, ThreadFeed, StreamPapers, TrendPapers, \
    ThreadFeedThreads


logger = logging.getLogger(__name__)

User = get_user_model()


@app.task(ignore_result=True)
def reset_all_main_streams():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        reset_stream.delay(user_pk)


@app.task(ignore_result=True)
def update_all_main_streams():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_stream.delay(user_pk)


@app.task(ignore_result=True)
def reset_all_main_trends():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        reset_trend.delay(user_pk)


@app.task(ignore_result=True)
def update_all_main_trends():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_trend.delay(user_pk)


@app.task()
def update_stream(user_pk, name='main'):
    """Async task / Update user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=name)
    # update
    feed.update()
    return user_pk


@app.task()
def reset_stream(user_pk, name='main'):
    """Async task / Reset user feed"""
    # create/update main feed
    feed, _ = Stream.objects.get_or_create(user_id=user_pk, name=name)
    # reset
    feed.clear_all()
    # update
    feed.update()
    return user_pk


@app.task()
def update_trend(user_pk, name='main'):
    trend, _ = Trend.objects.get_or_create(user_id=user_pk, name=name)
    trend.update()
    return user_pk


@app.task()
def reset_trend(user_pk, name='main'):
    trend, _ = Trend.objects.get_or_create(user_id=user_pk, name=name)
    trend.update()
    return user_pk


@app.task(ignore_result=True)
def reset_all_main_threadfeed():
    us_pk = User.objects.annotate(social_count=Count(F('social_auth')))\
        .exclude(social_count__lt=1)\
        .values_list('id', flat=True)
    for user_pk in us_pk:
        update_threadfeed.delay(user_pk)


@app.task()
def update_threadfeed(user_pk, name='main'):
    df, _ = ThreadFeed.objects.get_or_create(user_id=user_pk, name=name)
    df.update()
    return user_pk


@app.task(ignore_result=True)
def populate_stream(res, stream_id):
    """Populate stream based on res"""

    stream = Stream.objects.get(id=stream_id)

    logger.info('Populating stream {id}'.format(id=stream.id))
    stream.set_state('ING')

    # update score threshold
    stream.score_threshold = stream.user.settings.stream_score_threshold

    if res:     # res can be empty is user library is empty
        # reformat
        res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                        for r in res if r['score'] > stream.score_threshold])
        pids = list(res_dic.keys())

        # clean stream
        stream.clean_not_in(pids)

        # Update existing StreamPapers
        with transaction.atomic():
            sp_update = StreamPapers.objects.select_for_update()\
                .filter(stream=stream, paper_id__in=pids)
            update_pids = []
            for sp in sp_update:
                sp.score = res_dic[sp.paper_id]['score']
                sp.save()
                update_pids.append(sp.paper_id)

        # Create new StreamPapers
        with transaction.atomic():
            create_objs = []
            create_pids = set(pids).difference(set(update_pids))
            # Cross-check that thread id is still alive
            if create_pids:
                query = 'SELECT id FROM library_paper WHERE id IN (VALUES{0})'\
                    .format(','.join(['({0})'.format(x) for x in create_pids]))
                create_pids = [p.id for p in Paper.objects.raw(query)]

            for id_ in create_pids:
                create_objs.append(StreamPapers(
                    stream=stream,
                    paper_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            StreamPapers.objects.bulk_create(create_objs)

        # Get last user visit. use for the tagging matched with 'new'
        try:
            last_seen = LastSeen.objects.when(user=stream.user)
            StreamPapers.objects.filter(created__lt=last_seen).update(new=False)
        except LastSeen.DoesNotExist:
            pass

    stream.updated_at = timezone.now()
    stream.save(update_fields=('updated_at', ))
    stream.set_state('IDL')
    logger.info('Updating stream {id} done'.format(id=stream.id))


@app.task(ignore_result=True)
def populate_trend(res, trend_id):
    """Populate trend based on res"""
    # TODO: Factorize populate_trend and populate_stream

    trend = Trend.objects.get(id=trend_id)

    logger.info('Populating trend {id}'.format(id=trend.id))
    trend.set_state('ING')

    # update score threshold
    trend.score_threshold = trend.user.settings.stream_score_threshold

    if res:     # res can be empty is user library is empty
        # reformat
        res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                        for r in res if r['score'] > trend.score_threshold])
        pids = list(res_dic.keys())

        # clean trend
        trend.clean_not_in(pids)

        # Update existing TrendPapers
        with transaction.atomic():
            tp_update = TrendPapers.objects.select_for_update()\
                .filter(trend=trend, paper_id__in=pids)
            update_pids = []
            for tp in tp_update:
                tp.score = res_dic[tp.paper_id]['score']
                tp.save()
                update_pids.append(tp.paper_id)

        # Create new TrendPapers
        with transaction.atomic():
            create_objs = []
            create_pids = set(pids).difference(set(update_pids))

            # Cross-check that thread id is still alive
            if create_pids:
                query = 'SELECT id FROM library_paper WHERE id IN (VALUES{0})'\
                    .format(','.join(['({0})'.format(x) for x in create_pids]))
                create_pids = [p.id for p in Paper.objects.raw(query)]

            for id_ in create_pids:
                create_objs.append(TrendPapers(
                    trend=trend,
                    paper_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            TrendPapers.objects.bulk_create(create_objs)

        # Get last user visit. use for the tagging matched with 'new'
        try:
            last_seen = LastSeen.objects.when(user=trend.user)
            TrendPapers.objects.filter(created__lt=last_seen).update(new=False)
        except LastSeen.DoesNotExist:
            pass

    trend.updated_at = timezone.now()
    trend.save(update_fields=('updated_at', ))

    trend.set_state('IDL')
    logger.info('Updating trend {id} done'.format(id=trend.id))


@app.task(ignore_result=True)
def populate_threadfeed(res, threadfeed_id):
    """Populate threadfeed"""
    # TODO: Factorize populate_trend and populate_stream

    threadfeed = ThreadFeed.objects.get(id=threadfeed_id)

    logger.info('Updating thread feed {id}'.format(id=threadfeed.id))
    threadfeed.set_state('ING')

    # update score threshold
    threadfeed.score_threshold = threadfeed.user.settings.threadfeed_score_threshold

    if res:     # res can be empty is user library is empty
        # reformat
        res_dic = dict([(r['id'], {'score': r['score'], 'date': r['date']})
                        for r in res if r['score'] > threadfeed.score_threshold])
        tids = list(res_dic.keys())

        # clean threadfeed
        threadfeed.clean_not_in(tids)

        # Update existing TrendFeedThreads
        with transaction.atomic():
            tfp_update = ThreadFeedThreads.objects.select_for_update()\
                .filter(threadfeed=threadfeed, thread_id__in=tids)
            update_tids = []
            for tfp in tfp_update:
                tfp.score = res_dic[tfp.thread_id]['score']
                tfp.save()
                update_tids.append(tfp.thread_id)

        # Create new TrendFeedThreads
        create_objs = []
        create_pids = set(tids).difference(set(update_tids))
        with transaction.atomic():

            # Cross-check that thread id is still alive
            if create_pids:
                query = 'SELECT id FROM threads_thread WHERE id IN (VALUES{0})'\
                    .format(','.join(['({0})'.format(x) for x in create_pids]))
                create_pids = [t.id for t in Thread.objects.raw(query)]

            for id_ in create_pids:
                create_objs.append(ThreadFeedThreads(
                    threadfeed=threadfeed,
                    thread_id=id_,
                    score=res_dic[id_]['score'],
                    date=res_dic[id_]['date'],
                    new=True))
            ThreadFeedThreads.objects.bulk_create(create_objs)

        # Get last user visit. use for the tagging matched with 'new'
        try:
            last_seen = LastSeen.objects.when(user=threadfeed.user)
            ThreadFeedThreads.objects.filter(created__lt=last_seen).update(new=False)
        except LastSeen.DoesNotExist:
            pass

    threadfeed.updated_at = timezone.now()
    threadfeed.save(update_fields=('updated_at', ))

    threadfeed.set_state('IDL')
    logger.info('Updating thread feed {id} done'.format(id=threadfeed.id))
