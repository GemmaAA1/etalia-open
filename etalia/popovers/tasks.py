# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from config.celery import celery_app as app
from .models import UserPopOverUpdateDisplay, PopOver, UserPopOver

logger = logging.getLogger(__name__)


@app.task()
def update_popovers_display(user_pk):
    upoud = UserPopOverUpdateDisplay.objects.get(user_id=user_pk)
    upoud.update_display()


@app.task()
def init_popovers(user_pk):
    pos = PopOver.objects.all()
    for po in pos:
        UserPopOver.objects.get_or_create(user_id=user_pk, popover=po)
    upoud, _ = UserPopOverUpdateDisplay.objects.get_or_create(user_id=user_pk)
    upoud.update_display()
    return user_pk
