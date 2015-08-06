import re
import string
import random
from dateutil.parser import parse
import datetime

from ..constants import MENDELEY_PT, ZOTERO_PT
from core.parsers import Parser


class ParserBackend(Parser):
    """Backend Parser

    Attributes:
        user_info_template: A dictionary that stores data on how paper
            is related to user based on provider info

    """
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
        """Return a random ID"""
        return ''.join(random.choice(chars) for _ in range(size))


class ParserMendeley(ParserBackend):
    """Mendeley Parser"""

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        # Journal
        journal['title'] = entry.source

        if isinstance(entry.identifiers, dict):
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
        if isinstance(entry.identifiers, dict):
            for key, val in entry.identifiers.items():
                if key == 'doi':
                    paper['id_doi'] = val
                elif key == 'pmid':
                    paper['id_pmi'] = val
                elif key == 'pii':
                    paper['id_pii'] = val
                elif key == 'arxiv':
                    paper['id_arx'] = val
                    # overwrite paper type
                    paper['type'] = 'PRE'
                    paper['publish_status'] = 'preprint'
                else:
                    if not key == 'issn':
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

        try:
            user_info['created'] = parse(str(entry.created) or 'Nothing')
        except ValueError:
            pass
        try:
            user_info['last_modified'] = parse(str(entry.last_modified)
                                               or 'Nothing')
        except ValueError:
            pass
        user_info['starred'] = entry.starred
        user_info['authored'] = entry.authored

        return user_info


class ParserZotero(ParserBackend):
    """Zotero Parser"""

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        journal['title'] = entry.get('publicationTitle', '')

        issn_s = [a.strip() for a in entry.get('ISSN', '').split(',')]
        journal['id_issn'] = issn_s[0]
        if len(issn_s) > 1:
            journal['id_eissn'] = issn_s[1]

        # arxiv case
        journal['id_arx'] = self.process_arxiv(entry, 'journal')

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # Type
        # TODO UPDATE THIS
        type_ = entry.get('itemType')
        paper['type'] = dict(ZOTERO_PT).get(type_, '')

        # Title
        paper['title'] = entry.get('title', '')

        # Identifiers
        paper['id_doi'] = entry.get('DOI', '')

        # special case arxiv
        paper['id_arx'] = self.process_arxiv(entry, 'paper')
        if paper['type'] == 'BOO':
            paper['id_isbn'] = entry.get('ISBN', '')

        # special case book
        paper['id_arx'] = self.process_arxiv(entry, 'paper')
        if paper['id_arx']:
            paper['type'] = 'PRE'
            paper['publish_status'] = 'preprint'

        # generate id_oth is no doi
        if not any([paper[key] for key in ['id_doi', 'id_pmi', 'id_pii',
                                           'id_arx', 'id_oth', 'id_isbn']]):
            paper['id_oth'] = 'rand{0}'.format(self.id_oth_generator())

        # URL
        if entry.get('url', ''):
            paper['url'] = entry.get('url', '')
        elif paper['id_doi']:
            paper['url'] = '{base_url}{doi}'.format(
                base_url='http://doi.org/',
                doi=paper['id_doi'])
        else:
            paper['url'] = ''

        # Published date
        try:
            paper['date_pp'] = parse(entry.get('date',
                                               'Nothing') or 'Nothing').date()
        except ValueError:
            pass

        # Volume, issue, page
        paper['volume'] = entry.get('volume', '')
        paper['issue'] = entry.get('issue', '')
        paper['page'] = entry.get('pages', '')

        # Language

        # Abstract
        paper['abstract'] = entry.get('abstractNote', '')

        return paper

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.get('creators', '')
        for auth in full_authors:
            author = self.author_template.copy()
            author['first_name'] = auth.get('firstName', '')
            author['last_name'] = auth.get('lastName', '')
            authors.append(author)

        return authors

    def parse_corp_authors(self, entry):
        corp_authors = []
        return corp_authors

    def parse_user_info(self, entry):

        user_info = self.user_info_template.copy()
        try:
            user_info['created'] = parse(str(entry.get('dateAdded', ''))
                                         or 'Nothing')
        except ValueError:
            pass
        try:
            user_info['last_modified'] = parse(
                str(entry.get('dateLastmodified', '')) or 'Nothing')
        except ValueError:
            pass

        return user_info

    @staticmethod
    def process_arxiv(entry, return_what):
        """Retrieve arxiv specific fields"""

        paper_id_arx = ''
        journal_id_arx = ''
        if entry.get('libraryCatalog', '') == 'arXiv.org':
            str_proc = entry.get('publicationTitle', '')
            res = re.match(r'arXiv:(?P<id>[\d\.]+)\s\[(?P<sub>[\w-]+)\]',
                           str_proc)

            if not res: # try this one 'arXiv:cond-mat/9910006' (depreciated)
                res = re.match(r'arXiv:(?P<sub>[\w-]+)\/(?P<id>[\d\.]+)',
                               str_proc)

            if not res: # try that one 'arXiv:1206.1773' (depreciated)
                res = re.match(r'arXiv:(?P<id>[\d\.]+)',
                               str_proc)

            if res:
                paper_id_arx = res.groupdict().get('id', '')
                journal_id_arx = 'arxiv.'+res.groupdict().get('sub', '')
            elif entry.get('reportNumber', ''):
                paper_id_arx = entry.get('reportNumber', '')
                journal_id_arx = 'arxiv.common'
            elif entry.get('url', ''):
                res = re.match(r'http:\/\/arxiv.org\/abs\/(?P<id>[\d\.]+)',
                               entry.get('url', ''))
                paper_id_arx = res.groupdict().get('id', '')
                journal_id_arx = 'arxiv.abs'

        if return_what == 'paper':
            return paper_id_arx
        elif return_what == 'journal':
            return journal_id_arx
