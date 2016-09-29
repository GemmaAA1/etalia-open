# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import sitemaps
from django.core.urlresolvers import reverse


class StaticViewSitemap(sitemaps.Sitemap):
    priority = 0.5
    changefreq = 'weekly'
    protocol = 'https'

    def items(self):
        return ['core:home',
                'core:about',
                'core:contact',
                'core:help',
                'core:support',
                'core:terms_privacy',
                'core:terms_use']

    def location(self, item):
        return reverse(item)


class TrendSitemap(sitemaps.Sitemap):
    priority = 0.5
    changefreq = 'daily'
    protocol = 'https'

    def items(self):
        return ['library:papers']

    def location(self, item):
        return reverse(item)