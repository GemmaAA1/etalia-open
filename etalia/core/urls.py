# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf.urls import url
from . import views

from django.contrib.sitemaps.views import sitemap
from .sitemaps import StaticViewSitemap

sitemaps = {
    'static': StaticViewSitemap,
}

urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^about/$', views.about, name='about'),
    url(r'^terms-of-use/$', views.terms_use, name='terms_use'),
    url(r'^terms-of-privacy/$', views.terms_privacy, name='terms_privacy'),
    url(r'^support/$', views.support, name='support'),
    url(r'^press/$', views.press, name='press'),
    url(r'^press/(?P<slug>[a-zA-Z0-9-]+)_(?P<pk>[0-9]+)/$', views.press_slug, name='press-slug'),
    url(r'^press/(?P<pk>[0-9]+)/$', views.press_pk, name='press-pk'),
    url(r'^contact/$', views.contact, name='contact'),
    url(r'^help/$', views.help, name='help'),
    url(r'^test-failing-task$', views.test_failing_task),
    url(r'^sitemap\.xml$', sitemap, {'sitemaps': sitemaps},
        name='django.contrib.sitemaps.views.sitemap')
]
