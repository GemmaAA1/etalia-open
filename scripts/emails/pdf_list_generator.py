# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from etalia.library.models import Paper, Journal

# Arxiv
ps = Paper.objects.filter(id_arx__isnull=False)
urls = []
for p in ps:
    urls.append('http://arxiv.org/pdf/{id}v1.pdf'.format(id=p.id_arx))
file = open('/home/ubuntu/arxiv_pdf_list.txt', 'w+')
for item in urls:
    file.write('{0}\n'.format(item))
file.close()


# Plos One
j = Journal.objects.get(title__iexact='Plos One')
ps = Paper.objects.filter(journal=j)
urls = []
for p in ps:
    if p.id_doi:
        urls.append('http://journals.plos.org/plosone/article/asset?id={doi}.PDF'.format(doi=p.id_doi))
file = open('/home/ubuntu/plosone_pdf_list.txt', 'w+')
for item in urls:
    file.write('{0}\n'.format(item))
file.close()


#