"""To retrieve a list of affiliation from pubmed"""

import re
from Bio import Entrez
from Bio import Medline
from dateutil import parser


query = [ "abstract"]

Entrez.email = 'nicolas.pannetier@gmail.com'

# Fetch results
affs = []
ret_start = 0
ret_max = 200
while True:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=ret_max, restart=ret_start)
    record = Entrez.read(handle)
    handle.close()
    idlist = record["IdList"]
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",
                           retmode="text")
    records = Medline.parse(handle)
    affs += [reco.get('AD', '') for reco in records]
    ret_start += ret_max

    if len(affs > 10000):
        break


pre_affs = []
for record in affs:
    tmp = re.findall(r'([\w\s\-,]+)[\.|;]\s', record)
    pre_affs += [tm for tm in tmp if len(tm) > 20]
pre_affs = list(set(pre_affs))

