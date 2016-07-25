# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, patterns, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()
router.register(r'feeds/streampapers', views.StreamPapersViewSets)
router.register(r'feeds/trendpapers', views.TrendPapersViewSets)
router.register(r'feeds/threadfeedthreads', views.ThreadFeedThreadsViewSets)
router.register(r'feeds/stream', views.StreamViewSets)
router.register(r'feeds/trend', views.TrendViewSets)
router.register(r'feeds/threadfeed', views.ThreadFeedViewSets)

urlpatterns = [
    url(r'', include(router.urls)),
]
