import re
import abc
from library.models import Paper
from library.forms import PaperForm
from .constants import PUBMED_PT
from abc import ABCMeta, abstractmethod


class Converter(metaclass=ABCMeta):

    __metaclass__ = abc.ABCMeta

    item = {
        'type': '',
        'id_doi': '',
        'id_arx': '',
        'id_pii': '',
        'id_pmi': '',
        'id_oth': '',
        'title': '',
        'abstract': '',
        'journal': None,
        'volume': '',
        'issue': '',
        'page': '',
        'date_ep': None,
        'date_p': None,
        'date_lr': None,
        'url': '',
        'language': '',
        'is_aip': False,
        'is_pre_print': False,
    }

    @abc.abstractmethod
    def convert(self, entry):
        return NotImplementedError


class PubmedConverter(Converter):

    TEMPLATE_IDS = {'id_doi': r'(.+)\s\[doi\]',
                    'id_pii': r'(.+)\s\[pii\]'}

    def convert(self, entry):

        # Type
        try:
            type_ = entry.get('PT', '')
            self.item['type'] = dict(PUBMED_PT)[type_]
        except KeyError:
            self.item['type'] = ''

        # Identifiers
        # match template
        for id_ in entry.get('AID', ''):
            for key, pattern in self.TEMPLATE_IDS.items():
                match = re.match(pattern, id_)
                if match:
                    self.item[key] = match.groups()[0]
        # pmid
        self.item['id_pmi'] = entry.get('PMID', '')

        # Article in press
        if self.item.get('PST', '') == 'aheadofprint':
            self.item['is_aip'] = True

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