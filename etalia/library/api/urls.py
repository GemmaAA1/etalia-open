# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, patterns, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()
router.register(r'papers', views.PaperViewSet)
router.register(r'journals', views.JournalViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]