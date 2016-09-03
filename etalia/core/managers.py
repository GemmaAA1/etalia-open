# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import requests
from requests.exceptions import HTTPError, ConnectionError, Timeout
import feedparser
from habanero import Crossref
from habanero.exceptions import RequestError
from config.celery_settings.development import CELERY_ALWAYS_EAGER
from Bio import Entrez, Medline
from django.conf import settings
from django.db import IntegrityError, transaction
from django.db.models import Q
from django.contrib.auth import get_user_model
from etalia.library.models import Paper, Author, AuthorPaper, Journal, \
    CorpAuthor, CorpAuthorPaper, PaperUser
from etalia.users.models import UserLibPaper
from etalia.library.forms import PaperForm, AuthorForm, JournalForm, CorpAuthorForm
from etalia.feeds.models import StreamPapers, TrendPapers
from etalia.consumers.parsers import CrossRefPaperParser, ArxivPaperParser, \
    PubmedPaperParser, ElsevierPaperParser
from etalia.core.parsers import PaperParser
from etalia.threads.forms import PubPeerCommentForm, PubPeerForm, ThreadForm
from etalia.threads.models import Thread, PubPeer, PubPeerComment
import logging


logger = logging.getLogger(__name__)


def concatenate_last_names(authors):
    return ' '.join([a['last_name'] for a in authors])


class ObjDict(dict):
    """ ObjDict class provides a mean to treat dict as object with attributes"""

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class ConsolidateManager(object):
    """ The ConsolidateManager class provides methods for consolidate entries
    as return by parsers (e.g. ZoteroParser, MendeleyParser, etc)

    The consolidation workflow depends on the field found in the entry.
    The consolidation is based on Cross Ref, Pubmed and Arxiv databse.
    """

    METHOD_UPDATE_FIELDS = {
        'crossref': ['id_doi', 'volume', 'issue', 'page', 'url', 'title',
                     'date_ep', 'date_pp', 'type'],
        'pubmed': ['id_pmi', 'id_doi', 'volume', 'issue', 'page', 'url',
                   'title', 'date_ep', 'date_pp', 'abstract', 'type'],
        'arxiv': ['id_arx', 'volume', 'issue', 'page', 'url', 'title',
                  'date_ep', 'date_pp', 'date_lr', 'abstract', 'type'],
        'elsevier': ['id_doi', 'id_pii', 'volume', 'issue', 'page', 'url',
                     'title', 'date_ep', 'date_pp', 'abstract', 'type']
    }

    def __init__(self, entry):
        self.entry = entry

    def consolidate(self):
        """ConsolidateManager dispatcher workflow """
        self.entry['is_trusted'] = False
        if self.entry['paper'].get('id_arx'):
            self.consolidate_with('arxiv')
        elif self.entry['paper'].get('id_pmi'):
            self.consolidate_with('pubmed')
        else:
            seq = ['crossref', 'pubmed', 'elsevier', 'arxiv']
            count = 0
            while (self.entry.get('is_trusted') == False
                or self.entry.get('abstract') == '') and count < len(seq):
                self.consolidate_with(seq[count])
                count += 1

        return self.entry

    def update_entry(self, new, method_name):
        """Update entry"""
        fields = self.METHOD_UPDATE_FIELDS[method_name]
        self.update_entry_paper_fields(new, fields)
        self.entry['authors'] = new['authors']
        self.entry['journal'] = new['journal']
        self.entry['corp_authors'] = new['corp_authors']

    def update_entry_paper_fields(self, new, fields):
        """Update entry based on method fields allowed"""
        for k in fields:
            self.entry['paper'][k] = new['paper'][k]

    def check_query_match(self, new):
        """Check that new entry against entry attribute
        Check only title and last name of authors"""
        pattern = re.compile('[\W_]+')
        # compare title
        if not pattern.sub('', new['paper']['title'].lower()) == \
                pattern.sub('', self.entry['paper']['title'].lower()):
            return False

        # compare last name authors
        fetch_names = set([a['last_name'].lower() for a in new['authors']])
        store_names = set([a['last_name'].lower() for a in self.entry['authors']])
        if not store_names.issubset(fetch_names):
            return False

        return True

    def consolidate_with(self, method_name):
        """ConsolidateManager method dispatcher"""

        method = getattr(self, 'get_{name}'.format(name=method_name))
        if method_name == 'arxiv':
            doc_id = self.entry['paper'].get('id_arx')
        else:
            doc_id = self.entry['paper'].get('id_doi')
        query = requests.utils.quote('{title} {authors}'.format(
            title=self.entry['paper'].get('title', ''),
            authors=concatenate_last_names(self.entry.get('authors', []))))
        query = requests.utils.quote(query, safe='')
        new_entry = method(doc_id=doc_id, query=query)
        if new_entry and (doc_id or self.check_query_match(new_entry)):
            self.update_entry(new_entry, method_name)
            self.entry['id_gen'] = ''
            self.entry['is_trusted'] = True
        else:
            self.entry['is_trusted'] = False
        return self.entry

    @staticmethod
    def get_pubmed(doc_id='', query=''):
        """Return data from pubmed api"""
        logger.info('Starting new HTTP connection (1): pubmed')
        try:
            parser = PubmedPaperParser()
            email = settings.CONSUMER_PUBMED_EMAIL
            Entrez.email = email
            doi = doc_id
            if doi:
                handle = Entrez.esearch(db='pubmed',
                                        term='{doi}[AID]'.format(doi=doi))
            else:
                handle = Entrez.esearch(db='pubmed',
                                        term=query)

            record = Entrez.read(handle)
            handle.close()
            id_list = list(record["IdList"])
            handle = Entrez.efetch(db="pubmed",
                                   id=id_list,
                                   rettype="medline",
                                   retmode="text")
            records = Medline.parse(handle)
            entries = [record for record in records]
            if entries and entries[0].get('PMID'):
                entry = entries[0]
                return parser.parse(entry)
        except IOError:
            pass

        return None

    @staticmethod
    def get_crossref(doc_id='', query=''):
        """Return data from crossref api"""
        try:
            parser = CrossRefPaperParser()
            cr = Crossref()
            doi = doc_id
            if doi:
                try:
                    entry = cr.works(ids=[doi]).get('message')
                    return parser.parse(entry)
                except HTTPError:
                    pass
            entries = cr.works(query=query, limit=1).get('items')
            if entries:
                entry = entries[0]
                return parser.parse(entry)
        except:
            pass

        return None

    @staticmethod
    def get_arxiv(doc_id='', query=''):
        """Return data from arXiv api"""
        try:
            url_query = 'http://export.arxiv.org/api/query?search_query='

            if doc_id:
                id_arx = re.sub(r'v[0-9]+', '', doc_id)
                if id_arx.startswith('arXiv:'):
                    id_arx = id_arx.strip('arXiv:')
                resp = requests.get('{url}{id}'.format(url=url_query,
                                                       id=id_arx))
            else:
                resp = requests.get('{url}{query}&max_results=1'.format(
                    url=url_query,
                    query=query)
                )

            entries = feedparser.parse(resp.text).get('entries')
            if entries:
                parser = ArxivPaperParser()
                entry = entries[0]
                return parser.parse(entry)
        except Timeout or HTTPError or ConnectionError:
            pass

        return None

    @staticmethod
    def get_elsevier(doc_id='', query=''):
        api_key = settings.CONSUMER_ELSEVIER_API_KEY
        url_query = 'http://api.elsevier.com/content/search/index:SCIDIR?query='
        headers = {'X-ELS-APIKey': api_key}
        try:
            data = {
                'count': 1,
                'field': u'doi,coverDate,coverDisplayDate,url,identifier,title,'
                         'publicationName,issueIdentifier,coverDisplayName,'
                         'authors,creator,description,startingPage,issn,'
                         'endingPage',
                'start': 0,
            }
            if doc_id:
                q = '{base}DOI({id})'.format(base=url_query, id=doc_id)
            else:
                q = '{base}{query}'.format(base=url_query, query=query)
            resp = requests.post(q, data=data, headers=headers)
            if 'search-results' in resp.json().keys():
                entries = resp.json().get('search-results').get('entry')
                parser = ElsevierPaperParser()
                entry = entries[0]
                return parser.parse(entry)
        except Timeout or HTTPError or ConnectionError:
            pass
        return None


class PaperManager(object):

    def __init__(self, consolidate=False):
        self.consolidate = consolidate

    def get_or_create_from_entry(self, entry, fetch_journal=True):
        """Create or Retrieve paper based on entry data"""
        if self.entry_is_valid(entry):
            # Paper, (Journal if any)
            paper, journal = self.get_or_create_paper_from_entry(entry['paper'])
            if paper:
                if not paper.is_trusted:
                    # Journal
                    if fetch_journal and not journal:
                        paper.journal = \
                            self.get_journal_from_entry(entry['journal'])
                    # Authors
                    self.get_or_create_authors_from_entry(
                        entry['authors'],
                        entry['corp_authors'],
                        paper
                    )

                    paper.is_trusted = entry['is_trusted']
                    paper.save()

                    # consolidate async
                    # NB: This task if run as eager true conflict with
                    # update_lib pipeline because it modifies paper id.
                    if self.consolidate and not CELERY_ALWAYS_EAGER:
                        from etalia.consumers.tasks import consolidate_paper
                        consolidate_paper.apply_async(args=(paper.id, ),
                                                      countdown=30)

                    # Embed async
                    paper.embed()

            return paper, journal
        return None, None

    def consolidate_paper(self, paper, force=False):

        # convert paper to entry type
        entry = {
            'paper': PaperForm(instance=paper).initial,
            'authors': [AuthorForm(instance=auth).initial for auth in
                        paper.authors.all().order_by('authorpaper')],
            'journal': JournalForm(instance=paper.journal).initial,
            'corp_authors': [CorpAuthorForm(instance=ca).initial for ca in
                            paper.corp_author.all()]
        }
        # update None id_* to blank (for database query)
        for k, v in entry['paper'].items():
            if k.startswith('id_') and v is None:
                entry['paper'][k] = ''
        for k, v in entry['journal'].items():
            if k.startswith('id_') and v is None:
                entry['journal'][k] = ''

        # consolidate entry
        new_entry = ConsolidateManager(entry).consolidate()
        ep = new_entry['paper']

        # check for uniqueness
        # maybe the new_entry is already in DB, let's check
        qs = Paper.objects\
            .filter(
                Q(id_doi=ep['id_doi']) |
                Q(id_pmi=ep['id_pmi']) |
                Q(id_pii=ep['id_pii']) |
                Q(id_arx=ep['id_arx']) |
                Q(id_isbn=ep['id_isbn']) |
                Q(id_oth=ep['id_oth']))\
            .exclude(id=paper.id)

        if qs.count() > 0:
            # if so relink related objects to same_paper
            self.update_related_paper_fk(paper.id, qs[0].id)
            # delete current paper
            paper.delete()
            paper = qs[0]
        else:
            try:
                paper = qs[0]
            except IndexError:
                pass

        if force or not paper.is_trusted:
            # update paper with new_entry
            paper = self.update_paper_from_entry(new_entry, paper)
            paper.is_trusted = new_entry['is_trusted']
            paper.save()

            # send embedding
            paper.embed()

        return paper

    def update_paper_from_entry(self, entry, paper):
        ep = entry['paper']
        form = PaperForm(ep, instance=paper)
        if form.is_valid():
            paper = form.save()

        if entry['journal']:
            journal = self.get_journal_from_entry(entry['journal'])
            paper.journal = journal
            paper.save(update_fields=['journal'])

        # Clean authors and rewrite from trusted source
        paper.authors.all().delete()
        self.get_or_create_authors_from_entry(entry['authors'],
                                              entry['corp_authors'],
                                              paper)
        return paper

    @staticmethod
    def update_related_paper_fk(from_id, to_id):
        from etalia.threads.models import Thread
        # Try/except to deal with unique_together constraint
        try:
            UserLibPaper.objects.filter(paper_id=from_id).update(paper_id=to_id)
        except IntegrityError:
            pass
        try:
            PaperUser.objects.filter(paper_id=from_id).update(paper_id=to_id)
        except IntegrityError:
            pass
        try:
            StreamPapers.objects.filter(paper_id=from_id).update(paper_id=to_id)
        except IntegrityError:
            pass
        try:
            TrendPapers.objects.filter(paper_id=from_id).update(paper_id=to_id)
        except IntegrityError:
            pass
        try:
            Thread.objects.filter(paper_id=from_id).update(paper_id=to_id)
        except IntegrityError:
            pass

    @staticmethod
    def get_or_create_paper_from_entry(ep):
        """Create or Retrieve paper"""
        try:
            paper = Paper.objects.get(
                Q(id_doi=ep['id_doi']) |
                Q(id_pmi=ep['id_pmi']) |
                Q(id_pii=ep['id_pii']) |
                Q(id_arx=ep['id_arx']) |
                Q(id_isbn=ep['id_isbn']) |
                Q(id_oth=ep['id_oth'])
            )
            if paper.is_trusted:
                return paper, paper.journal
        except Paper.MultipleObjectsReturned:
            paper = Paper.objects.filter(
                Q(id_doi=ep['id_doi']) |
                Q(id_pmi=ep['id_pmi']) |
                Q(id_pii=ep['id_pii']) |
                Q(id_arx=ep['id_arx']) |
                Q(id_isbn=ep['id_isbn']) |
                Q(id_oth=ep['id_oth'])
            ).first()

        except Paper.DoesNotExist:
            paper = None

        form = PaperForm(ep, instance=paper)
        if form.is_valid():
            paper = form.save()
            return paper, None
        return None, None

    def get_journal_from_entry(self, ej):
        """Retrieve Journal. Either with ID or exact title"""
        if self.journal_has_id(ej):
            try:
                # NB: 1st and 2nd conditions are because most of providers
                # do not distinguish e-issn and issn
                journal = Journal.objects.get(
                    Q(id_issn=ej['id_issn']) |
                    Q(id_eissn=ej['id_issn']) |
                    Q(id_eissn=ej['id_eissn']) |
                    Q(id_arx=ej['id_arx']) |
                    Q(id_oth=ej['id_oth']))
            except Journal.DoesNotExist:
                journal = None
            except Journal.MultipleObjectsReturned:
                # clean up journals
                journals = Journal.objects.filter(
                    Q(id_issn=ej['id_issn']) |
                    Q(id_eissn=ej['id_issn']) |
                    Q(id_eissn=ej['id_eissn']) |
                    Q(id_arx=ej['id_arx']) |
                    Q(id_oth=ej['id_oth']))
                # keep journals with most reference and delete other
                counts = [(j, j.paper_set.count()) for j in journals]
                counts = sorted(counts, key=lambda x: x[1])
                journal, _ = counts.pop()
                for i in counts:
                    Paper.objects.filter(journal_id=i[0].id)\
                        .update(journal_id=journal.id)
                    Journal.objects.get(id=i[0].id).delete()
                pass
        else:
            try:
                journal = Journal.objects.get(
                    title__iexact=ej['title'])
            except Journal.DoesNotExist:
                journal = None
                pass
            except Journal.MultipleObjectsReturned:
                # clean up journals
                journal = Journal.objects.filter(
                    Q(id_issn=ej['id_issn']) |
                    Q(id_eissn=ej['id_issn']) |
                    Q(id_eissn=ej['id_eissn']) |
                    Q(id_arx=ej['id_arx']) |
                    Q(id_oth=ej['id_oth'])).first()
                pass
        return journal

    @staticmethod
    def get_or_create_authors_from_entry(ea, eca, paper):
        """Create Author relationship"""
        with transaction.atomic():
            # Get or Create authors
            for pos, item_author in enumerate(ea):
                author, _ = Author.objects.get_or_create(
                    first_name=item_author['first_name'],
                    last_name=item_author['last_name'])
                AuthorPaper.objects.get_or_create(
                    paper=paper,
                    author=author,
                    position=pos)

            # create/get corp author
            for pos, item_corp_author in enumerate(eca):
                corp_author, _ = CorpAuthor.objects.get_or_create(
                    name=item_corp_author['name']
                )
                CorpAuthorPaper.objects.get_or_create(
                    paper=paper,
                    corp_author=corp_author)

    @staticmethod
    def entry_is_valid(entry):
        if entry['paper'].get('title', '') \
                and entry['authors'] \
                and entry['paper'].get('type', ''):
            return True
        return False

    @staticmethod
    def journal_has_id(item_journal):
        return any([True for key, val in item_journal.items()
                    if key.startswith('id_') and val])


class PubPeerManager(object):
    """Pipe PubPeer data to etalia"""

    def get_or_create_related_paper(self, doi):
        try:
            paper = Paper.objects.get(id_doi=doi)
        except Paper.DoesNotExist:
            paper_template = PaperParser.paper_template.copy()
            paper_template['id_doi'] = doi
            entry = {'paper': paper_template}
            new_entry = ConsolidateManager(entry).consolidate()
            new_entry['paper']['source'] = 'PPR'
            paper, _ = PaperManager().get_or_create_from_entry(new_entry)
        return paper

    def add_or_update_entry(self, entry):

        User = get_user_model()
        thread_entry = entry['thread']
        pubpeer_entry = entry['pubpeer']
        comments_entry = entry['comments']

        # Thread
        try:
            thread = Thread.objects.get(pubpeer__pubpeer_id=
                                        pubpeer_entry['pubpeer_id'])
            thread_entry['paper'] = thread.paper.id
            thread_entry['title'] = thread.title
            thread_entry['user'] = thread.user.id
        except Thread.DoesNotExist:
            thread = None
            doi = pubpeer_entry['doi']
            paper = self.get_or_create_related_paper(doi)
            if paper:
                thread_entry['paper'] = paper.id
                thread_entry['title'] = 'Comment on: {0}'.format(paper.title)
                thread_entry['user'] = User.objects.get(
                    email=settings.CONSUMER_PUBPEER_USER_EMAIL
                ).id
            else:
                return None

        form = ThreadForm(thread_entry, instance=thread)
        if form.is_valid():
            thread = form.save()
            thread.published_at = thread_entry['published_at']
            thread.save(update_fields=['published_at'])
        else:
            raise ValueError('ThreadForm is invalid {0}'
                             .format(form._errors))

        # PubPeer
        if thread:
            try:
                pb = PubPeer.objects.get(
                    pubpeer_id=pubpeer_entry['pubpeer_id']
                )
                pubpeer_entry['thread'] = pb.thread_id
            except PubPeer.DoesNotExist:
                pb = None
                pubpeer_entry['thread'] = thread.id
            form = PubPeerForm(pubpeer_entry, instance=pb)
            if form.is_valid():
                pb = form.save()
            else:
                raise ValueError('PubPeerForm is invalid {0}'
                                 .format(form._errors))

            # PubPeerComment
            with transaction.atomic():
                for c in comments_entry:
                    try:
                        pbc = PubPeerComment.objects.get(
                            pubpeercomment_id=c['pubpeercomment_id']
                        )
                        c['pubpeer'] = pbc.pubpeer_id
                    except PubPeerComment.DoesNotExist:
                        pbc = None
                        c['pubpeer'] = pb.id
                    form = PubPeerCommentForm(c, instance=pbc)
                    if form.is_valid():
                        form.save()
                    else:
                        raise ValueError('PubPeerComment form is invalid {0}'
                                         .format(form._errors))
        thread.pubpeer.update_is_active()
        return thread
