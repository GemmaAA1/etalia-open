# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import re
import json

input_file = '/Users/nicolaspannetier/Projects/paperstream/docs/J_Entrez.txt'
# input_file = '/Users/nicolaspannetier/Projects/paperstream/docs/test.txt'
output_file = '/Users/nicolaspannetier/Projects/paperstream/docs/J_Entrez_parse.json'

template = {
    'title':        '',
    'short_title':  '',   # Iso abbreviated title when available
    'id_issn':      '',   # Print identifier
    'id_eissn':     '',   # Online/Electronic identifier
    'id_arx':       '',   # Arxiv identifier
    'id_oth':       '',   # Other identifier
    'period':       '',   # publishing period recurrence
    'language':     '',   # language of the journal
    'url':          '',   # url of the journal page
    'publisher':    '',   # publisher name
    'scope':        ''    #
}

# regular expression to find
regexp = {
    'title': 'JournalTitle: (?P<short_title>[\w\s\'\-]+)\n',
    'short_title': 'IsoAbbr: (?P<short_title>[\w\s\'\-]+)\n',
    'id_issn': 'ISSN \(Print\): (?P<id_issn>[\w\d\-]+)\n',
    'id_eissn': 'ISSN \(Online\): (?P<id_issn>[\w\d\-]+)\n',
}

# Process
N = 8
data = []
count = 0
with open(input_file, 'r') as infile:
    lines = infile.readlines()
    num_blocks = len(lines) // N
    for i in range(num_blocks):
        block_str = ''.join([line for line in lines[i*N:(i+1)*N]])
        tmp = template.copy()
        for key, reg in regexp.items():
            res = re.findall(reg, block_str)
            if res:
                tmp[key] = res[0]
        data.append(tmp)

        if not i % 100:
            print('# of entry processed: {0}'.format(i))

# Dump to json output_file
json.dump(data, open(output_file, 'w'), indent=2, sort_keys=True)




