import abc
from .models import Paper, Journal, Author, CorpAuthor
from .forms import PaperForm, JournalForm, AuthorForm, CorpAuthorForm
from abc import ABCMeta, abstractmethod


class Parser(metaclass=ABCMeta):

    __metaclass__ = abc.ABCMeta

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
        journal = self.parse_journal(entry)
        paper = self.parse_paper(entry)
        authors = self.parse_authors(entry)
        corp_authors = self.parse_corp_authors(entry)

        return {'journal': journal, 'authors': authors, 'paper': paper,
                'corp_authors': corp_authors}