#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


def warm_up():
    pe_dispatcher.delay('dummy')
    te_dispatcher.delay('dummy')
    nlp_dispatcher.delay('dummy')


if __name__ == '__main__':

    # import django
    # django.setup()
    from etalia.nlp.tasks import pe_dispatcher, te_dispatcher, nlp_dispatcher

    warm_up()
