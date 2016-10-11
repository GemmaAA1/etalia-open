# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
from django.utils import timezone
import requests
from time import sleep
from bs4 import BeautifulSoup


class SimpleCrawlerListMixin(object):

    DOMAIN = ''
    LIST_PATTERN = ''
    DETAIL_PATTERN = {'name': '', 'class': ''}
    TIMEOUT = 1

    def get_urls_list(self, page):
        # build list url and parse details link
        url = self.LIST_PATTERN.format(domain=self.DOMAIN, page=page)
        sleep(self.TIMEOUT)
        doc = requests.get(url).text
        soup = BeautifulSoup(doc, 'html.parser')

        # Get list of paper detail
        tas = soup.find_all(self.DETAIL_PATTERN.get('name'),
                            {'class': self.DETAIL_PATTERN.get('class')})
        links = []
        for ta in tas:
            href = ta['href']
            if not href.startswith(self.DOMAIN):
                href = self.DOMAIN + href
            links.append(href)

        return links

    def crawl_links(self, links):
        details = []
        for url in links:
            details.append(requests.get(url).text)
            sleep(self.TIMEOUT)
        return details

    def crawl_to_date(self, date):
        page = 1
        links = []
        # Crawl listing pages
        while not self.stop_crawl_date(links, date):
            links += self.get_urls_list(page)
            page += 1
        # Crawl detail page
        return self.crawl_links(links)

    def stop_crawl_date(self, links, date):
        if links:
            last_link = links[-1]
            res = re.match(r'.+/content\/early\/([\d]+)\/([\d]+)\/([\d]+).+', last_link)
            if res:
                current_date = timezone.datetime(int(res.groups()[0]),
                                                 int(res.groups()[1]),
                                                 int(res.groups()[2])).date()
                if current_date < date.date():
                    return True
            else:
                raise ValueError('regexp failed: {url}'.format(url=last_link))
        return False

