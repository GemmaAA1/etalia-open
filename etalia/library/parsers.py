# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from etalia.core.parsers import Parser
import datetime


class CrossRefParser(Parser):
    """CrossRef Parser
    """

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
        if len(j_titles) > 1:
            journal['short_title'] = j_titles[1]
        journal['title'] = j_titles[0]

        issns = entry.get('ISSN')
        if len(issns) > 1:
            journal['id_eissn'] = issns[1]
        journal['id_issn'] = issns[0]

        return journal

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # type
        type_ = entry.get('type', 'UNKNOWN')
        paper['type'] = dict(self.CROSSREF_PT).get(type_)

        # title
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
            paper['date_pp'] = datetime.date(y, m, d)

        # Volume, issue, page
        paper['volume'] = entry.get('volume')
        paper['issue'] = entry.get('issue')
        paper['page'] = entry.get('page')

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
        return self.corp_author_template.copy()
