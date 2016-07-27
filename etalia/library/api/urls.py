# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, patterns, url
from django.views.decorators.cache import cache_page, cache_control, never_cache

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()
router.register(r'library/papers', views.PaperViewSet)
router.register(r'library/journals', views.JournalViewSet)
router.register(r'library/authors', views.AuthorViewSet)
router.register(r'library/states', views.PaperStateViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]


