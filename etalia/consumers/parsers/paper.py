# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import datetime
from bs4 import BeautifulSoup
from dateutil.parser import parse
from nameparser import HumanName

from etalia.core.parsers import PaperParser


class PubmedPaperParser(PaperParser):
    """Pubmed PaperParser"""

    TEMPLATE_IDS = {'id_doi': r'(.+)\s\[doi\]',
                    'id_pii': r'(.+)\s\[pii\]'}

    PUBMED_PT = (
        ('JOURNAL ARTICLE', 'JOU'),
        ('LETTER',          'LET'),
        ('EDITORIAL',       'EDI'),
        ('NEWS',            'NEW'),
        ('CONGRESSES',      'PRO'),
        ('REVIEW',          'REV'),
        ('PATENTS',         'PAT'),
        ('UNKNOWN',         ''),
    )
    TYPE = 'PUB'

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        # Journal
        IS = entry.get('IS', '')
        issn_match = re.findall(r'(?P<issn1>\d{4}-\d{3}[\dX])\s\(Printing\)|'
                                r'(?P<issn2>\d{4}-\d{3}[\dX])\s\(Linking\)|'
                                r'(?P<eissn>\d{4}-\d{3}[\dX])\s\(Electronic\)',
                                IS)
        # the following two lines sucks because the regexp matching above sucks
        issn = ''
        eissn = ''
        if issn_match:
            try:
                issn = [i[0] or i[1] for i in issn_match if i[0] or i[1]][0]
                eissn = [i[2] for i in issn_match if i[2]][0]
            except IndexError:
                pass

        journal['id_issn'] = issn
        journal['id_eissn'] = eissn
        journal['title'] = entry.get('JT', '')

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # Type
        type_ = entry.get('PT', [''])
        if type_:
            types = [dict(self.PUBMED_PT).get(typ.upper(), '') for typ in type_]
            types = [t for t in types if t]
            if types:
                paper['type'] = types[0]

        # Identifiers
        # match template
        ids = entry.get('AID', [''])
        if ids:
            for id_ in ids:
                for key, pattern in self.TEMPLATE_IDS.items():
                    match = re.match(pattern, id_)
                    if match:
                        paper[key] = match.groups()[0]
        # pmid
        paper['id_pmi'] = entry.get('PMID', '')

        # URL
        if paper['id_doi']:
            paper['url'] = '{base_url}{doi}'.format(
                base_url='http://doi.org/',
                doi=paper['id_doi'])
        elif paper['id_pmi']:
            paper['url'] = '{base_url}{pmid}'.format(
                base_url='http://http://www.ncbi.nlm.nih.gov/pubmed/',
                pmid=paper['id_pmi'])
        else:
            paper['url'] = ''

        # Paper published status
        paper['publish_status'] = entry.get('PST', '')

        # Published date
        try:
            paper['date_ep'] = parse(entry.get('DEP', 'a_date_that_excepts'))
        except ValueError or AttributeError:
            paper['date_ep'] = None
        try:
            paper['date_pp'] = parse(entry.get('DP', 'a_date_that_excepts'))
        except ValueError or AttributeError:
            paper['date_pp'] = None
        try:
            paper['date_lr'] = parse(entry.get('LR', 'a_date_that_excepts'))
        except ValueError or AttributeError:
            paper['date_lr'] = None

        # Volume, issue, page
        paper['volume'] = entry.get('VI', '')
        paper['issue'] = entry.get('IP', '')
        paper['page'] = entry.get('PG', '')

        # Language
        lang = entry.get('LA', [''])
        if lang:
            paper['language'] = lang[0].upper()

        # Abstract
        paper['abstract'] = entry.get('AB', '')

        # Title
        paper['title'] = entry.get('TI', '')

        return paper

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.get('FAU', [''])
        for auth in full_authors:
            author = self.author_template.copy()
            name = HumanName(auth)
            author['last_name'] = name.last.strip()
            author['first_name'] = name.first.strip()
            if name.middle.strip():
                author['first_name'] += ' ' + name.middle.strip()
            authors.append(author)
        return authors

    def parse_corp_authors(self, entry):

        corp_authors = []

        corp_authors_items = entry.get('CN', '')
        if corp_authors_items:
            if isinstance(corp_authors_items, list):
                for corp_auth in corp_authors_items:
                    corp_author = self.corp_author_template.copy()
                    corp_author['name'] = corp_auth
                    corp_authors.append(corp_author)
            else:
                corp_author = self.corp_author_template.copy()
                corp_author['name'] = corp_authors_items
                corp_authors.append(corp_author)
        return corp_authors


class ArxivPaperParser(PaperParser):
    """Arxiv PaperParser"""

    TYPE = 'ARX'

    def parse_authors(self, entry):
        authors = []

        full_authors = entry.get('authors', [''])
        for auth in full_authors:
            author = self.author_template.copy()
            name = HumanName(auth['name'])
            author['last_name'] = name.last.strip()
            author['first_name'] = name.first.strip()
            if name.middle:
                author['first_name'] += ' ' + name.middle.strip()
            authors.append(author)

        return authors

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        prim_cat = entry.get('arxiv_primary_category', '')
        if prim_cat:
            term = entry.get('arxiv_primary_category').get('term', '')
            if term:
                journal['id_arx'] = 'arxiv.{term}'.format(term=term)

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        paper['type'] = 'PRE'
        paper['publish_status'] = 'preprint'

        patterns = [
            # For (old?) article already published?
            'http://arxiv.org/abs/(?P<jid>[\w-]+)/((?P<pid>[\d\.]+)(v\d)?)',
            # Current pattern
            'http://arxiv.org/abs/((?P<pid>[\d\.]+)(v\d)?)'
        ]
        for pattern in patterns:
            res = re.search(pattern, entry.get('id', ''))
            if res:
                paper['id_arx'] = '/'.join(res.groupdict().values())
                break

        paper['url'] = entry.get('link', '')

        try:
            paper['date_ep'] = parse(entry.get('published',
                                               'a_date_that_excepts'))
        except ValueError or AttributeError:
            paper['date_ep'] = None
        try:
            paper['date_lr'] = parse(entry.get('updated',
                                               'a_date_that_excepts'))
        except ValueError or AttributeError:
            paper['date_lr'] = None
        paper['abstract'] = ' '.join(entry.get('summary', '').split())

        paper['title'] = ' '.join(entry.get('title', '').split())

        return paper

    def parse_corp_authors(self, entry):
        return []


class ElsevierPaperParser(PaperParser):
    """Elsevier PaperParser"""

    TYPE = 'ELS'

    def parse_authors(self, entry):
        authors = []

        auths = entry.get('authors')
        if auths:
            full_authors = auths.get('author')
            for auth in full_authors:
                author = self.author_template.copy()
                try:
                    author['first_name'] = auth.get('given-name', '').strip()
                    author['last_name'] = auth.get('surname', '').strip()
                    authors.append(author)
                except AttributeError:
                    pass

        return authors

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        issn = entry.get('prism:issn', '')
        if issn:
            id_issn = '{0}-{1}'.format(issn[:4], issn[4:])
            journal['id_issn'] = id_issn

        journal['title'] = entry.get('prism:publicationName', '')

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # type
        paper['type'] = 'JOU'

        if 'Available online' in entry.get('prism:coverDisplayDate', ''):
            paper['publish_status'] = 'aheadofprint'
            date_ep = re.sub(
                r'Available online',
                '',
                entry.get('prism:coverDisplayDate', 'rejected')).strip()
            try:
                paper['date_ep'] = parse(date_ep)
            except ValueError or AttributeError:
                paper['date_ep'] = None
        else:
            paper['publish_status'] = 'ppublish'
            try:
                paper['date_pp'] = parse(entry.get('prism:coverDisplayDate',
                                                   'a_date_that_excepts'))
            except ValueError or AttributeError:
                paper['date_pp'] = None
            paper['page'] = '{start}-{end}'.format(
                start=entry.get('prism:startingPage', ''),
                end=entry.get('prism:endingPage', ''))
            paper['volume'] = entry.get('prism:volume', '')
            paper['issue'] = entry.get('prism:issueIdentifier', '')

        paper['id_doi'] = entry.get('prism:doi', '')

        # pii
        if entry.get('prism:url'):
            paper['id_pii'] = entry.get('prism:url').split('/')[-1]

        paper['issue'] = entry.get('prism:issueIdentifier', '')

        paper['url'] = '{urlbase}{pii}'.format(
            urlbase='http://www.sciencedirect.com/science/article/pii/',
            pii=paper['id_pii'])

        paper['abstract'] = entry.get('dc:description', '')

        paper['title'] = entry.get('dc:title', '')

        return paper

    def parse_corp_authors(self, entry):
        return []


class CrossRefPaperParser(PaperParser):
    """CrossRef PaperParser
    """

    TYPE = 'CRO'

    CROSSREF_PT = (
        ('journal-article',         'JOU'),
        ('report',                  'JOU'),
        ('book',                    'BOO'),
        ('book-set',                'BOO'),
        ('book-series',             'BOO'),
        ('editor-book',             'BOO'),
        ('reference-book',          'BOO'),
        ('book-track',              'BOs'),
        ('book-part',               'BOS'),
        ('proceedings-article',     'PRO'),
        ('proceedings',             'PRO'),
        ('patent',                  'PAT'),
        ('dissertation',            'THE'),
        ('other',                   ''),
        ('UNKNOWN',                 ''),
    )

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        # Journal
        j_titles = entry.get('container-title', [])
        if j_titles:
            if len(j_titles) > 1:
                journal['short_title'] = j_titles[1]
            try:
                journal['title'] = j_titles[0]
            except IndexError:
                journal['title'] = ''
        else:
            # specific case for bioRxiv based on doi match registrant
            doi = entry.get('DOI')
            if re.match(r'10.1101/\d+', doi):
                journal['id_oth'] = 'biorxiv'

        issns = entry.get('ISSN')
        if issns:
            if len(issns) > 1:
                journal['id_eissn'] = issns[1]
            try:
                journal['id_issn'] = issns[0]
            except IndexError:
                journal['id_issn'] = ''

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # type
        type_ = entry.get('type', 'UNKNOWN')
        paper['type'] = dict(self.CROSSREF_PT).get(type_)

        # title
        if entry.get('title'):
            paper['title'] = entry.get('title')[0]

        # publisher
        publisher = entry.get('publisher')

        # id
        paper['id_doi'] = entry.get('DOI')
        if publisher.lower().startswith('elsevier'):
            paper['id_pii'] = entry.get('alternative-id', [''])[0]

        # url
        paper['url'] = entry.get('URL')

        # Published date
        if 'published-print' in entry:
            date = entry.get('published-print')['date-parts'][0]
            y = date[0]
            m = date[1] if len(date) > 1 else 1
            d = date[2] if len(date) > 2 else 1
            paper['date_pp'] = datetime.date(y, m, d)
        if 'published-online' in entry:
            date = entry.get('published-online')['date-parts'][0]
            y = date[0]
            m = date[1] if len(date) > 1 else 1
            d = date[2] if len(date) > 2 else 1
            paper['date_ep'] = datetime.date(y, m, d)

        # Volume, issue, page
        paper['volume'] = entry.get('volume', '')
        paper['issue'] = entry.get('issue', '')
        paper['page'] = entry.get('page', '')

        return paper

    def parse_authors(self, entry):
        authors = []
        for auth in entry.get('author', []):
            author = self.author_template.copy()
            author['last_name'] = auth.get('family')
            author['first_name'] = auth.get('given', '')
            authors.append(author)
        return authors

    def parse_corp_authors(self, entry):
        return []


class BiorxivPaperParser(PaperParser):
    """Parse HTML detail page of a paper as found for ex at
    http://biorxiv.org/content/early/2016/10/04/050245
    """

    TYPE = 'BIO'

    def pre_process(self, entry):
        return BeautifulSoup(entry, 'html.parser')

    def parse_journal(self, entry):

        journal = self.journal_template.copy()
        journal['short_title'] = 'BioRxiv'
        journal['title'] = 'BioRxiv'
        journal['id_oth'] = 'biorxiv'

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        paper['type'] = 'PRE'
        if entry.h1:
            paper['title'] = entry.h1.text

        # id
        doi_block = entry.find('span', {'class': 'highwire-cite-metadata-doi'})
        if doi_block:
            res = re.match(r'.+http://dx.doi.org/([\.\w\d\/]+)', doi_block.text)
            if res:
                paper['id_doi'] = res.groups()[0]

        # Published date
        date_block = entry.find('li', {'class': 'published'})
        if date_block:
            res = re.match(r'Posted ([\d\w\s,]+)', date_block.text)
            if res:
                paper['date_ep'] = parse(res.groups()[0])

        # Abstract
        abstract_block = entry.find('p', {'id': 'p-2'})
        if abstract_block:
            paper['abstract'] = abstract_block.text

        # Url
        if paper['date_ep'] and paper['id_doi']:
            paper['url'] = 'http://biorxiv.org/content/early/{y}/{m}/{d}/{id}'.format(
                y=paper['date_ep'].year,
                m=paper['date_ep'].month,
                d=paper['date_ep'].day,
                id=paper['id_doi'].split('/')[-1]
            )
        elif paper['id_doi']:
            paper['url'] = 'http://dx.doi.org/{doi}'.format(doi=paper['id_doi'])

        return paper

    def parse_authors(self, entry):

        author_block = entry\
            .find('div', {'class': 'pane-highwire-article-citation'})\
            .findAll('span', {'class': 'highwire-citation-author', })

        authors = []
        for auth in author_block:
            author = self.author_template.copy()
            ln_block = auth.find('span', {'class': 'nlm-given-names'})
            if ln_block:
                author['last_name'] = ln_block.text
            fn_block = auth.find('span', {'class': 'nlm-surname'})
            if fn_block:
                author['first_name'] = fn_block.text
            authors.append(author)
        return authors

    def parse_corp_authors(self, entry):
        return []
