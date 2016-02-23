# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
import csv

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

# input_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_thomsonreuters_cleaned_light.csv'
# output_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_thomsonreuters_cleaned_light.json'
# input_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_thomsonreuters_cleaned.csv'
# output_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_thomsonreuters_cleaned.json'
# input_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_arxiv_cleaned.csv'
# output_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_arxiv_cleaned.json'
input_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_medline_cleaned.csv'
output_file = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/static/populate/journals/20150510_medline_cleaned.json'

# read csv file
with open(input_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\n')
    lines = []
    for row in reader:
        lines.append(row[0])

# process
data = []
keys = lines[0].split(';')
for j, line in enumerate(lines[1:]):
    tmp = template.copy()
    row_data = line.split(';')
    for i, key in enumerate(keys):
        tmp[key] = row_data[i]
    data.append(tmp)

    if not j % 100:
        print('# of entry processed: {0}'.format(j))

# Dump to json output_file
json.dump(data, open(output_file, 'w'), indent=2, sort_keys=True)

