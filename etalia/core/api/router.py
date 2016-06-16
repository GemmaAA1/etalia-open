# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework.routers import SimpleRouter, DefaultRouter
from django.conf import settings


class SharedAPIRootRouter(SimpleRouter):
    shared_router = DefaultRouter()

    def register(self, *args, **kwargs):
        self.shared_router.register(*args, **kwargs)
        super().register(*args, **kwargs)


def api_urls():
    """ Scan installed_apps and register api urls"""
    from importlib import import_module
    for app in settings.INSTALLED_APPS:
        try:
            import_module(app + '.api.urls')
        except (ImportError, AttributeError):
            pass
    return SharedAPIRootRouter.shared_router.urls