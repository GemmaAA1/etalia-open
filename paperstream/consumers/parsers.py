import re
from dateutil.parser import parse

from .constants import PUBMED_PT
from library.parsers import Parser


class ParserPubmed(Parser):

    TEMPLATE_IDS = {'id_doi': r'(.+)\s\[doi\]',
                    'id_pii': r'(.+)\s\[pii\]'}

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
            paper['type'] = \
                [dict(PUBMED_PT).get(typ.upper(), '') for typ in type_][0]

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


class ParserArxiv(Parser):

    def parse_authors(self, entry):
        authors = []

        full_authors = entry.get('authors', [''])
        for auth in full_authors:
            author = self.author_template.copy()
            # assuming all first name are only initialized [e.g E. J. Martin]
            split_auth = list(map(str.strip, auth['name'].split('.')))
            if len(split_auth) > 1:
                first = ' '.join(split_auth[:-1])
                last = split_auth[-1]
            else:
                last = auth['name'].strip()
                first = ''
            author['last_name'] = last
            author['first_name'] = first
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

        paper['id_arx'] = re.sub(r'http://arxiv.org/abs/','', entry.get('id', ''))

        paper['url'] = entry.get('link', '')

        paper['date_ep'] = parse(entry.get('published', 'a_date_that_excepts'))
        paper['date_lr'] = parse(entry.get('updated', 'a_date_that_excepts'))

        paper['abstract'] = entry.get('summary', '')

        paper['title'] = entry.get('title', '')

        return paper

    def parse_corp_authors(self, entry):

        corp_author = self.corp_author_template.copy()

        return []


class ParserElsevier(Parser):

    def parse_authors(self, entry):
        authors = []

        full_authors = entry.get('authors', {'': ''}).get('author', [''])
        for auth in full_authors:
            author = self.author_template.copy()
            try:
                author['first_name'] = auth.get('given-name', '')
                author['last_name'] = auth.get('surname', '')
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
            date_ep = re.sub(r'Available online',
                             '',
                             entry.get('prism:coverDisplayDate', 'rejected')).strip()
            paper['date_ep'] = parse(date_ep)
        else:
            paper['publish_status'] = 'ppublish'
            paper['date_pp'] = parse(entry.get('prism:coverDisplayDate',
                                               'rejected'))
            paper['page'] = '{start}-{end}'.format(
                start=entry.get('prism:startingPage', ''),
                end=entry.get('prism:endingPage', ''))
            paper['volume'] = entry.get('prism:volume', '')
            paper['issue'] = entry.get('prism:issueIdentifier', '')

        paper['id_doi'] = entry.get('prism:doi', '')

        # pii
        paper['id_pii'] = re.sub('http://api.elsevier.com/content/article/pii:',
                                 '',
                                 entry.get('prism:url', ''))

        paper['issue'] = entry.get('prism:issueIdentifier', '')

        paper['url'] = '{urlbase}{pii}X'.format(
            urlbase='http://www.sciencedirect.com/science/article/pii/',
            pii=paper['id_pii'])

        paper['abstract'] = entry.get('dc:description', '')

        paper['title'] = entry.get('dc:title', '')


        return paper

    def parse_corp_authors(self, entry):

        corp_author = self.corp_author_template.copy()

        return []