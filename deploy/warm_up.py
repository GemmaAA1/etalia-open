#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


def warm_up():
    # We are looping here to deal with worker concurrency.
    # I.e warm up tasks are broadcasting to workers but
    # within one worker if concurrency is >1, task is only seen by 1 thread
    for _ in range(10):
        warmup_paper_engine.delay()
        warmup_thread_engine.delay()
        warmup_nlp_model.delay()


if __name__ == '__main__':

    # import django
    # django.setup()
    from etalia.nlp.tasks import warmup_paper_engine, warmup_thread_engine, \
        warmup_nlp_model

    warm_up()
