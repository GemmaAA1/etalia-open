# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import requests
from bs4 import BeautifulSoup
from .parsers import BiorxivPaperParser

domain = 'http://biorxiv.org/'
list_url = '{domain}/archive?field_highwire_a_epubdate_value[value]&page={page}'

page = 1
doc = requests.get(list_url.format(domain=domain, page=page)).text
soup = BeautifulSoup(doc, 'html.parser')

# Get paper links in list
tas = soup.findAll('a', {'class': 'highwire-cite-linked-title'})
links = []
for ta in tas:
    links.append(ta['href'])

parser = BiorxivPaperParser()
for d in links:
    doc = requests.get('{domain}{d}'.format(domain=domain, d=d)).text
    paper = parser.parse(doc)





