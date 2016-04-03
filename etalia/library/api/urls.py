# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, patterns, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()
router.register(r'library/papers', views.PaperViewSet)
router.register(r'library/journals', views.JournalViewSet)
router.register(r'library/authors', views.AuthorViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]