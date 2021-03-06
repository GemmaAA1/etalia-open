# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import string
import random
from nameparser import HumanName
from dateutil.parser import parse
import datetime
from etalia.core.parsers import PaperParser


class PaperParserBackend(PaperParser):
    """Backend PaperParser

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
        out = super(PaperParserBackend, self).parse(entry)
        return dict({'user_info': self.parse_user_info(entry)}, **out)

    @staticmethod
    def id_oth_generator(size=10, chars=string.ascii_uppercase + string.digits):
        """Return a random ID"""
        return ''.join(random.choice(chars) for _ in range(size))


class ParserOrcid(PaperParser):

    ORCID_PT = (
        ('JOURNAL_ARTICLE',         'JOU'),
        ('REPORT',                  'JOU'),
        ('BOOK',                    'BOO'),
        ('BOOK_CHAPTER',            'BOS'),
        ('RESEARCH_TOOL',           'JOU'),
        ('CONFERENCE_ABSTRACT',     'PRO'),
        ('CONFERENCE_PAPER',        'PRO'),
        ('PATENT',                  'PAT'),
        ('DISSERTATION',            'THE'),
        ('UNKNOWN',                 ''),
        ('',                        ''),
    )

    TYPE = 'ORC'

    def parse_authors(self, entry):

        authors = []

        contributors = entry.get('work-contributors')
        if contributors:
            auths = contributors.get('contributor')
            if auths:
                for auth in auths:
                    author = self.author_template.copy()
                    try:
                        name = auth.get('credit-name', {}).get('value', '')
                        hn = HumanName(name)
                        author['first_name'] = (hn.first + ' ' + hn.middle).strip()
                        author['last_name'] = hn.last
                        authors.append(author)
                    except AttributeError:
                        pass

        return authors

    def parse_journal(self, entry):
        journal = self.journal_template.copy()

        ids = entry.get('work-external-identifiers').get('work-external-identifier')
        for id_ in ids:
            if id_.get('work-external-identifier-type') == 'ISSN':
                journal['id_issn'] = id_.get('work-external-identifier-id', {}).get('value', '')

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        paper['type'] = dict(self.ORCID_PT).get(entry.get('work-type', ''))

        ids = entry.get('work-external-identifiers').get('work-external-identifier')
        for id_ in ids:
            if id_.get('work-external-identifier-type') == 'DOI':
                paper['id_doi'] = id_.get('work-external-identifier-id', {}).get('value', '')

        paper['title'] = entry.get('work-title', {}).get('title', {}).get('value', '')

        try:
            publication_date = entry.get('publication-date')
            year, month, day = None, 1, 1
            if publication_date:
                if publication_date.get('year'):
                    year = int(publication_date.get('year').get('value'))
                if publication_date.get('month'):
                    month = int(publication_date.get('month').get('value'))
                if publication_date.get('day'):
                    day = int(publication_date.get('day').get('value'))
                if year:
                    paper['date_pp'] = datetime.date(year, month, day)
        except TypeError:
            pass

        return paper

    def parse_corp_authors(self, entry):
        corp_authors = []
        return corp_authors


class ParserMendeley(PaperParserBackend):
    """Mendeley PaperParser"""

    MENDELEY_PT = (
        ('journal',                 'JOU'),
        ('book',                    'BOO'),
        ('book_section',            'BOS'),
        ('conference_proceedings',  'PRO'),
        ('working_paper',           'DRA'),
        ('patent',                  'PAT'),
        ('thesis',                  'THE'),
        ('UNKNOWN',                 ''),
    )

    TYPE = 'MEN'

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
        paper['type'] = dict(self.MENDELEY_PT).get(type_, '')

        # Title
        paper['title'] = entry.title.strip()

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
                        paper['id_oth'] = '{0}_{1}'.format(
                            key,
                            val[len(val) + len(key)-63:])

        if not any([paper[key] for key in ['id_doi', 'id_pmi', 'id_pii',
                                           'id_arx', 'id_oth']]):
            paper['id_oth'] = 'to-be-generated'

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
            if paper['id_arx']:  # this is an early view for sure
                paper['date_ep'] = datetime.date(year, month, day)
            else:
                paper['date_pp'] = datetime.date(year, month, day)

        # Volume, issue, page
        paper['volume'] = (entry.volume or '').strip()
        paper['issue'] = (entry.issue or '').strip()
        paper['page'] = (entry.pages or '').strip()

        # Language

        # Abstract
        paper['abstract'] = entry.abstract

        return paper

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.authors
        if full_authors:
            for auth in full_authors:
                author = self.author_template.copy()
                author['last_name'] = (auth.last_name or '').strip()
                author['first_name'] = (auth.first_name or '').strip()
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


class ParserZotero(PaperParserBackend):
    """Zotero PaperParser"""

    TYPE = 'ZOT'

    ZOTERO_PT = (
        ('journalArticle',  'JOU'),
        ('book',            'BOO'),
        ('bookSection',     'BOS'),
        ('conferencePaper', 'PRO'),
        ('thesis',          'THE'),
        ('patent',          'PAT'),
        ('letter',          'LET'),
        ('UNKNOWN',         ''),
    )

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
        type_ = entry.get('itemType')
        paper['type'] = dict(self.ZOTERO_PT).get(type_, '')

        # Title
        paper['title'] = entry.get('title', '').strip()

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
            paper['id_oth'] = 'to-be-generated'

        # URL
        if 'url' in entry:
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
        paper['volume'] = entry.get('volume', '').strip()
        paper['issue'] = entry.get('issue', '').strip()
        paper['page'] = entry.get('pages', '').strip()

        # Language

        # Abstract
        paper['abstract'] = entry.get('abstractNote', '').strip()

        return paper

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.get('creators', '')
        for auth in full_authors:
            author = self.author_template.copy()
            if auth.get('lastName', None):
                author['first_name'] = auth.get('firstName', '').strip()
                author['last_name'] = auth.get('lastName', '').strip()
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
            elif 'reportNumber' in entry:
                paper_id_arx = entry.get('reportNumber', '')
                journal_id_arx = 'arxiv.common'
            elif 'url' in entry:
                res = re.match(r'http:\/\/arxiv.org\/abs\/(?P<id>[\d\.]+)',
                               entry.get('url', ''))
                if res:
                    paper_id_arx = res.groupdict().get('id', '')
                    journal_id_arx = 'arxiv.abs'
                else:
                    res = re.match(r'http:\/\/arxiv.org\/abs\/(?P<jid>[\w-]+)\/(?P<id>[\d\.]+)',
                               entry.get('url', ''))
                    if res:
                        paper_id_arx = res.groupdict().get('id', '')
                        journal_id_arx = res.groupdict().get('jid', '')


        if return_what == 'paper':
            return paper_id_arx
        elif return_what == 'journal':
            return journal_id_arx
