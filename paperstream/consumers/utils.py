import re



def pubmed_to_pap(entry):

    item = {}

    # Identifiers
    # template for identifiers
    template_id = {'id_doi': r'(.+)\s\[doi\]',
                   'id_pii': r'(.+)\s\[pii\]'}
    # match
    for id_ in entry.get('AID', ''):
        for key, pattern in template_id.items():
            match = re.match(pattern, id_)
            if match:
                item[key] = match.groups()[0]
    # PMID
    item['id_pmi'] = entry.get('PMID', '')

    # Article in press
    if item.get('PST', '') == 'aheadofprint':
        item['is_aip'] = True
    else:
        item['is_aip'] = False

    # Published date


    # if entry.get('AID', '') and \
    #         entry.get('FAU', '') and \
    #         entry.get('AB', ''):



        # published date
        date = parse(entry.get('DP', '')).date() or \
               parse(entry.get('DA', '')).date()
        item['dat'] = date

        # url
        item['url'] = 'http://doi.org' + doi

        # issn
        rm = re.findall(r'(\d{4}-\d{3}[\dX])\s\(Electronic\)',
                        entry.get('IS', ''))
        if rm:
            item['e_issn'] = rm[0]
        rm = re.findall(r'(\d{4}-\d{3}[\dX])\s\(Linking\)',
                        entry.get('IS', ''))
        if rm:
            item['issn'] = rm[0]
        rm = re.match(r'(\d{4}-\d{3}[\dX])\s\(Printing\)',
                      entry.get('IS', ''))
        if rm:
            item['issn'] = rm[0]

        # journal title
        item['jou'] = entry.get('JT', '')

        # Authors
        item['aut'] = entry.get('FAU', '')

        # title
        item['tit'] = entry.get('TI', '')

        # abstract
        item['abs'] = entry.get('AB', '')

        items.append(item)

    np = 0
    if items:
    # re-format authors
        for item in items:
            tmp = item['aut'].copy()
            for i, auth in enumerate(tmp):
                item['aut'][i] = {}
                last, first = tuple(list(map(str.strip, tmp[i].split(','))))
                item['aut'][i]['last_name'] = last
                item['aut'][i]['first_name'] = first

        # Add new paper to database and count
        for item in items:
            paper, new = libu.add_paper(item)
            if new:
                np += 1