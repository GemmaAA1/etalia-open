#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


def warm_up():
    # We are looping here to deal with worker concurrency.
    # within one worker if concurrency is >1, task is only seen by 1 thread
    for _ in range(10):
        nlp_dispatcher.delay('dummy')
        pe_dispatcher.delay('dummy')
        te_dispatcher.delay('dummy')


if __name__ == '__main__':

    # import django
    # django.setup()
    from etalia.nlp.tasks import nlp_dispatcher, te_dispatcher, pe_dispatcher

    warm_up()
