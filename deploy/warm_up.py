#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from kombu.common import Broadcast


def warm_up():
    warmup_paper_engine.delay('dummy')
    warmup_thread_engine.delay('dummy')
    warmup_nlp_model.delay('dummy')


if __name__ == '__main__':

    # import django
    # django.setup()
    from etalia.nlp.tasks import warmup_paper_engine, warmup_thread_engine, \
        warmup_nlp_model

    warm_up()
