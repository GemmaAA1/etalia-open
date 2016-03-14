# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import abc
from etalia.library.models import Paper, Journal, Author, CorpAuthor
from etalia.library.forms import PaperForm, JournalForm, AuthorForm, \
    CorpAuthorForm
from abc import abstractmethod
import hashlib


class Parser(object):
    """Abstract Parser class

    Parser is used to parse entry in journal, paper, authors and corp_author
    so that it library database can be populated. This abstract class initialize
    several attribute based on their Model fields.

    Attributes:
        paper_template: Template for Paper
        journal_template: Template for Journal
        author_template: Template for Author
        corp_author_tempalte: Template for CorpAuthor
    """

    paper_template = dict([(field, Paper._meta.get_field(field).default)
                           for field in PaperForm.Meta.fields])

    journal_template = dict([(field, Journal._meta.get_field(field).default)
                             for field in JournalForm.Meta.fields])

    author_template = dict([(field, Author._meta.get_field(field).default)
                           for field in AuthorForm.Meta.fields])

    corp_author_template = \
        dict([(field, CorpAuthor._meta.get_field(field).default)
             for field in CorpAuthorForm.Meta.fields])

    @abstractmethod
    def parse_journal(self, entry):
        return NotImplementedError

    @abstractmethod
    def parse_paper(self, entry):
        return NotImplementedError

    @abstractmethod
    def parse_authors(self, entry):
        return NotImplementedError

    @abstractmethod
    def parse_corp_authors(self, entry):
        return NotImplementedError

    def parse(self, entry):
        """Parse entry

        Args:
            entry (dict): A dictionary of keys/values

        Returns:
            (dict): A dictionary with keys: 'journal', 'authors', 'paper',
                'corp_authors' and values as defined by class attributes
                (*_template)
        """
        journal = self.parse_journal(entry)
        paper = self.parse_paper(entry)
        authors = self.parse_authors(entry)
        corp_authors = self.parse_corp_authors(entry)

        # create id_oth based on parsed data if tagged accordingly
        if paper['id_oth'] == 'to-be-generated':
            paper['id_oth'] = 'gen' + self.generated_hash_id(paper,
                                                             journal,
                                                             authors)

        return {'journal': journal, 'authors': authors, 'paper': paper,
                'corp_authors': corp_authors}

    @staticmethod
    def generated_hash_id(paper, journal, authors):
        str_ = ''
        str_ += paper.get('title', '')
        str_ += ' '.join([author['last_name'] for author in authors])
        if 'title' not in journal:
            str_ += ''
        else:
            str_ += journal.get('title')
        return hashlib.sha1(str_.encode('utf-8')).hexdigest()