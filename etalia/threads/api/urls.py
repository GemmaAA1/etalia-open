# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, patterns, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()
router.register(r'thread/threads', views.ThreadViewSet)
router.register(r'thread/posts', views.ThreadPostViewSet)
router.register(r'thread/comments', views.ThreadCommentViewSet)
router.register(r'thread/states', views.ThreadUserViewSet)
router.register(r'thread/invites', views.ThreadUserInviteViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]