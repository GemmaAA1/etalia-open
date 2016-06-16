# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import include, url

from etalia.core.api.router import SharedAPIRootRouter
from . import views

router = SharedAPIRootRouter()

router.register(r'user/users', views.UserViewSet)
router.register(r'user/affiliations', views.AffiliationViewSet)
router.register(r'user/user-libs', views.UserLibViewSet)
router.register(r'user/user-lib-papers', views.UserLibPaperViewSet)
router.register(r'user/relationships', views.RelationshipViewSet)

urlpatterns = [
    url(r'', include(router.urls)),
]