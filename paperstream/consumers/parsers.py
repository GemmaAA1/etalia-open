import re
import abc
from library.models import Paper, Journal, Author, CorpAuthor
from library.forms import PaperForm, JournalForm, AuthorForm, CorpAuthorForm
from .constants import PUBMED_PT
from abc import ABCMeta, abstractmethod
from dateutil.parser import parse


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
    def parse_authors(self, entry):
        return NotImplementedError

    @abstractmethod
    def parse_journal(self, entry):
        return NotImplementedError

    @abstractmethod
    def parse_paper(self, entry):
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


class PubmedParser(Parser):

    TEMPLATE_IDS = {'id_doi': r'(.+)\s\[doi\]',
                    'id_pii': r'(.+)\s\[pii\]'}

    def parse_authors(self, entry):

        authors = []

        full_authors = entry.get('FAU', '')
        for auth in full_authors:
            author = self.author_template.copy()
            try:
                last, first = tuple(list(map(str.strip, auth.split(','))))
            except ValueError:
                last = auth.strip()
                first = ''
                pass
            author['last_name'] = last
            author['first_name'] = first
            authors.append(author)
        return authors

    def parse_corp_authors(self, entry):

        corp_authors = []

        corp_authors_items = entry.get('CN', '')
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

    def parse_paper(self, entry):

        paper = self.paper_template.copy()

        # Type
        try:
            type_ = entry.get('PT', '')
            paper['type'] = [dict(PUBMED_PT).get(typ, '') for typ in type_][0]
        except KeyError:
            paper['type'] = ''

        # Identifiers
        # match template
        for id_ in entry.get('AID', ''):
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
        else:
            paper['url'] = '{base_url}{pmid}'.format(
                base_url='http://http://www.ncbi.nlm.nih.gov/pubmed/',
                pmid=paper['id_pmi'])

        # Paper published status
        paper['publish_status'] = entry.get('PST', '')

        # Published date
        try:
            paper['date_ep'] = parse(entry.get('DEP', ''))
        except ValueError or AttributeError:
            paper['date_ep'] = None
        try:
            paper['date_pp'] = parse(entry.get('DP', ''))
        except ValueError or AttributeError:
            paper['date_ep'] = None
        try:
            paper['date_lr'] = parse(entry.get('LR', ''))
        except ValueError or AttributeError:
            paper['date_ep'] = None

        # Volume, issue, page
        paper['volume'] = entry.get('VI', '')
        paper['issue'] = entry.get('IP', '')
        paper['page'] = entry.get('PG', '')

        # Language
        paper['language'] = entry.get('LA', '')[0].upper()

        # Abstract
        paper['abstract'] = entry.get('AB', '')

        # Title
        paper['title'] = entry.get('TI')

        return paper

    def parse_journal(self, entry):

        journal = self.journal_template.copy()

        # Journal
        IS = entry.get('IS', '')
        issn_match = re.findall(r'(?P<issn1>\d{4}-\d{3}[\dX])\s\(Printing\)|'
                                r'(?P<issn2>\d{4}-\d{3}[\dX])\s\(Linking\)|'
                                r'(?P<eissn>\d{4}-\d{3}[\dX])\s\(Electronic\)',
                                IS)
        # the following two lines sucks because the regexp matching above sucks
        issn = [i[0] or i[1] for i in issn_match if i[0] or i[1]][0]
        eissn = [i[2] for i in issn_match if i[2]][0]
        # Use

        journal['id_issn'] = issn
        journal['id_eissn'] = issn
        journal['title'] = entry.get('JT', '')

        return journal