import re
import string
import random
from dateutil.parser import parse
import datetime

from ..constants import MENDELEY_PT
from library.parsers import Parser


class ParserBackend(Parser):

    user_info_template = {
        'created': None,
        'last_modified': None,
        'starred': None,
        'authored': None,
    }

    def parse_user_info(self, entry):
        raise NotImplementedError('Implemented in subclass')

    def parse(self, entry):
        out = super(ParserBackend, self).parse(entry)
        return dict({'user_info': self.parse_user_info(entry)}, **out)

    @staticmethod
    def id_oth_generator(size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))


class ParserMendeley(ParserBackend):

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        # Journal
        journal['title'] = entry.source

        for key, val in entry.identifiers.items():
            if key == 'issn':
                journal['id_issn'] = val
            if key == 'arxiv':
                journal['id_arx'] = 'arxiv.common'

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # Type
        type_ = entry.type
        paper['type'] = dict(MENDELEY_PT).get(type_, '')

        # Title
        paper['title'] = entry.title

        # Identifiers
        # match template
        ids = entry.identifiers
        for key, val in ids.items():
            if key == 'doi':
                paper['id_doi'] = val
            elif key == 'pmi':
                paper['id_pmi'] = val
            elif key == 'pii':
                paper['id_pii'] = val
            elif key == 'arxiv':
                paper['id_arx'] = val
                # overwrite paper type
                paper['type'] = 'PRE'
                paper['publish_status'] = 'preprint'
            else:
                paper['id_oth'] = '{0}_{1}'.format(key, val)

        if not any([paper[key] for key in ['id_doi', 'id_pmi', 'id_pii',
                                           'id_arx', 'id_oth']]):
            paper['id_oth'] = 'rand{0}'.format(self.id_oth_generator())

        # URL
        if entry.websites:
            tmp = entry.websites
            if isinstance(tmp, list):
                paper['url'] = entry.websites[0]
            else:
                paper['url'] = entry.websites
        elif paper['id_doi']:
            paper['url'] = '{base_url}{doi}'.format(
                base_url='http://doi.org/',
                doi=paper['id_doi'])
        elif paper['id_pmi']:
            paper['url'] = '{base_url}{pmid}'.format(
                base_url='http://http://www.ncbi.nlm.nih.gov/pubmed/',
                pmid=paper['id_pmi'])
        else:
            paper['url'] = ''

        # Published date
        day = entry.day or 1
        month = entry.month or 1
        year = entry.year

        if year:
            if paper['id_arx']:  # well this is an early view for sure
                paper['date_ep'] = datetime.date(year, month, day)
            else:
                paper['date_pp'] = datetime.date(year, month, day)

        # Volume, issue, page
        paper['volume'] = entry.volume
        paper['issue'] = entry.issue
        paper['page'] = entry.pages

        # Language

        # Abstract
        paper['abstract'] = entry.abstract

        return paper

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.authors
        for auth in full_authors:
            author = self.author_template.copy()
            author['last_name'] = auth.last_name
            author['first_name'] = auth.first_name
            authors.append(author)
        return authors

    def parse_corp_authors(self, entry):
        corp_authors = []
        return corp_authors

    def parse_user_info(self, entry):

        user_info = self.user_info_template.copy()

        user_info['created'] = parse(str(entry.created))
        user_info['last_modified'] = parse(str(entry.last_modified))
        user_info['starred'] = entry.starred
        user_info['authored'] = entry.authored

        return user_info


class ParserZotero(ParserBackend):
    pass