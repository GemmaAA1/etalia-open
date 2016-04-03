# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()

router.register(r'user/users', views.UserViewSet)
router.register(r'user/libraries', views.UserLibViewSet)
router.register(r'user/userlibpapers', views.UserLibPaperViewSet)
router.register(r'user/relationships', views.RelationshipViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]