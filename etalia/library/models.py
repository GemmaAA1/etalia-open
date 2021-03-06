# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.db import models, connection, transaction
from django.core.urlresolvers import reverse
from django.utils import timezone
from django.utils.text import slugify
from django.conf import settings
from django.contrib.sites.models import Site
from model_utils.fields import MonitorField

from etalia.core.models import TimeStampedModel, NullableCharField
from etalia.core.mixins import ModelDiffMixin

from .validators import validate_issn, validate_author_names
from .constants import LANGUAGES, PUBLISH_PERIODS, PAPER_TYPE, PUBLISH_STATUS, \
    PAPER_BANNED, PAPER_PINNED, PAPER_WATCH, PAPER_STORE, PAPER_ADDED, \
    PAPER_TRASHED
from .utils import langcode_to_langpap

from langdetect import detect


class Publisher(TimeStampedModel):
    """Publisher group
    """
    # Group name
    name = models.CharField(max_length=200, unique=True)
    # base url
    url = models.URLField(blank=True, default=True)

    def __str__(self):
        return self.name


class Journal(TimeStampedModel):
    """Periodicals
    """
    # NullableCharField is used to enforced uniqueness at the database level,
    # by default django save None charfield as '' which is considered as a string
    # in query and conflict with uniqueness
    id_issn = NullableCharField(max_length=9, blank=True, null=True,
                                default='', validators=[validate_issn],
                                unique=True, verbose_name='ISSN', db_index=True)
    id_eissn = NullableCharField(max_length=9, blank=True, null=True,
                                 default='', validators=[validate_issn],
                                 unique=True, verbose_name='e-ISSN',
                                 db_index=True)
    id_arx = NullableCharField(max_length=32, blank=True, null=True,
                               default='', unique=True, db_index=True,
                               verbose_name='Arxiv ID')
    id_oth = NullableCharField(max_length=32, blank=True, null=True,
                               default='', unique=True, db_index=True,
                               verbose_name='Other ID')

    # periodical title
    title = models.CharField(max_length=200, default='')
    # short title
    short_title = models.CharField(max_length=100, blank=True, default='')

    # publisher group
    publisher = models.ForeignKey(Publisher, null=True, default=None, blank=True)
    # url
    url = models.URLField(blank=True, default='')
    # Scope
    scope = models.TextField(blank=True, max_length=1000, default='')

    # language
    language = models.CharField(max_length=3, choices=LANGUAGES,
                                default='ENG', blank=True)

    # period
    period = models.CharField(max_length=200, choices=PUBLISH_PERIODS,
                              default='', blank=True)

    # is flag
    is_trusted = models.BooleanField(default=False)

    # is known from fixture
    is_in_fixture = models.BooleanField(default=True)

    # number of papers in journal
    lib_size = models.IntegerField(default=0)

    class Meta:
        ordering = ['title']

    def __str__(self):
        return self.print_short_title

    def get_absolute_url(self):
        return reverse('library:journal-slug',
                       kwargs={'pk': self.pk,
                               'slug': slugify(self.title)})

    def count_papers(self):
        self.lib_size = self.paper_set.count()
        self.save()
        return self.lib_size

    def count_papers_trusted(self):
        return self.paper_set.filter(is_trusted=True).count()

    def count_ids(self):
        count = 0
        # Loop through field starting with 'id_'
        for field in self._meta.fields:
            if field.name.startswith('id_') and getattr(self, field.name):
                count += 1
        return count

    @property
    def print_ids(self):
        ids_str = ''
        # Loop through field starting with 'id_'
        for field in self._meta.fields:
            if field.name.startswith('id_') and getattr(self, field.name):
                if ids_str:
                    ids_str += ', '
                ids_str += '{id}: {value}'.format(
                    id=field.verbose_name,
                    value=getattr(self, field.name))
        return ids_str

    @property
    def print_short_title(self):
        if len(self.short_title) > 2:
            return self.short_title
        else:
            if len(self.title) > 30:
                return self.title[:30]+'...'
            else:
                return self.title[:30]

    def print_title(self, trunc):
        if len(self.title) > trunc:
            return self.title[:]+'...'
        else:
            return self.title


class Author(TimeStampedModel):
    """Creators (authors) of matches
    """
    # TODO: add affiliation of authors ?

    # last name (capitalized)
    last_name = models.CharField(max_length=100,
                                 validators=[validate_author_names],
                                 default='')
    # first name (capitalize). first_name can store middle name
    first_name = models.CharField(max_length=100, blank=True, default='')
    # email
    email = models.EmailField(max_length=254, blank=True, default='')

    class Meta:
        unique_together = ('first_name', 'last_name')

    def __str__(self):
        if self.first_name:
            return self.first_name + ' ' + self.last_name
        else:
            return self.last_name

    @property
    def print_compact(self):
        # get initials
        if self.first_name:
            initials = '.'.join([name[0] for name in
                                 self.first_name.split(' ')]) + '.'
            return self.last_name + ' ' + initials
        else:
            return self.last_name

    @property
    def print_full(self):
        if self.first_name:
            return self.first_name + ' ' + self.last_name
        else:
            return self.last_name


class CorpAuthor(TimeStampedModel):

    name = models.CharField(max_length=128, blank=True, default='', unique=True)


class Paper(TimeStampedModel):
    """Scientific matches
    """

    # Published status
    publish_status = models.CharField(choices=PUBLISH_STATUS, blank=True,
                                      default='', max_length=20)
    publish_status_changed = MonitorField(monitor='publish_status')

    # Type of paper
    type = models.CharField(max_length=3, choices=PAPER_TYPE, blank=True,
                            default='')

    # identifiers (uniqueness defined thereafter)
    # NullableCharField is used to enforced uniqueness at the database level,
    # by default django save None charfield as '' which is considered as a string
    # in query and conflict with uniqueness
    id_doi = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='DOI',
                               db_index=True)
    id_arx = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='Arxiv',
                               db_index=True)
    id_pmi = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='PMID',
                               db_index=True)
    # not unique because publisher dependent
    id_pii = NullableCharField(max_length=64, blank=True, default='',
                               null=True, verbose_name='PII', db_index=True)
    id_isbn = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='ISBN',
                               db_index=True)
    id_oth = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='Other ID',
                               db_index=True)
    # Title
    title = models.CharField(max_length=500, blank=False, default='',
                             db_index=True)

    # Authors
    # authors
    authors = models.ManyToManyField(Author, through='AuthorPaper',
                                     blank=True, default=None, null=False)
    # corporate authorship
    corp_author = models.ManyToManyField(CorpAuthor, through='CorpAuthorPaper',
                                         blank=True, default=None)
    # Abstract
    abstract = models.TextField(blank=True, default='')

    # Journal
    journal = models.ForeignKey(Journal, null=True, blank=True, default=None,
                                db_index=True)
    # volume
    volume = models.CharField(max_length=200, blank=True, default='')
    # issue
    issue = models.CharField(max_length=200, blank=True, default='')
    # page
    page = models.CharField(max_length=200, blank=True, default='')

    # Key Terms
    key_terms = models.TextField(blank=True, default='')

    # Dates
    # date electronically published
    date_ep = models.DateField(null=True, blank=True, default=None,
                               db_index=True)
    # date of publication in issue journal
    date_pp = models.DateField(null=True, blank=True, default=None,
                               db_index=True)
    # date of paper last revised (e.g. arxiv, or publisher with e.g grant#)
    date_lr = models.DateField(null=True, blank=True, default=None)

    # date first seen
    date_fs = models.DateField(null=True, blank=True, db_index=True)

    # Consolidate date
    date_co = models.DateField(null=True, blank=True, db_index=True)

    # url where found
    url = models.URLField(blank=True, default='')
    # language
    language = models.CharField(max_length=3, choices=LANGUAGES,
                                default='ENG', blank=True)

    # Source
    source = models.CharField(blank=True, default='', max_length=20)
    source_changed = MonitorField(monitor='source')

    # Boolean
    # locked if all field have been verified or paper from 'trusted' source
    is_trusted = models.BooleanField(default=False, db_index=True)

    def save(self, **kwargs):
        if not self.date_fs:
            self.date_fs = timezone.now().date()
        if not self.date_co:
            dates = [self.date_ep, self.date_fs, self.date_pp]
            if any(dates):
                self.date_co = min([date for date in dates if date is not None])
        super(Paper, self).save(**kwargs)

    def get_absolute_url(self):
        return reverse('library:paper-slug', kwargs={'pk': self.pk,
                                                     'slug': slugify(self.title)})

    @property
    def short_title(self):
        if len(self.title) > 49:
            return '{0}...'.format(self.title[:50])
        else:
            return self.title

    @property
    def date(self):
        if self.date_ep:
            return self.date_ep
        elif self.date_pp:
            return self.date_pp
        else:
            return self.date_fs

    @property
    def print_compact(self):
        terms = [self.print_compact_authors,
                 self.title,
                 self.journal.print_short_title if self.journal else '',
                 str(self.date.year),
                 self.print_volume_issue_page]
        terms = [t for t in terms if t]
        return '. '.join(terms)

    @property
    def print_volume_issue_page(self):
        vol = str(self.volume)
        iss = str(self.issue)
        page = self.page
        if iss and vol and page:
            return '{vol}({iss}):{page}'.format(vol=vol, iss=iss, page=page)
        if self.issue and self.volume:
            return '{vol}({iss})'.format(vol=vol, iss=iss)
        if self.volume and self.page:
            return '{vol}:{page}'.format(vol=vol, page=page)
        if self.page:
            return self.page
        return ''

    @property
    def print_authors(self):
        return self.print_compact_authors()

    @property
    def print_compact_authors(self):
        authors = self.authors.all()
        authors_str = ''
        if authors:
            for author in authors:
                if authors_str:
                    authors_str += ', '
                authors_str += author.print_compact
            return authors_str
        else:
            return 'Unknown Authors'

    @property
    def print_full_authors(self):
        authors = self.authors.all()
        authors_str = ''
        if authors:
            for author in authors:
                if authors_str:
                    authors_str += ', '
                authors_str += author.print_full
            return authors_str
        else:
            return 'Unknown Authors'

    @property
    def print_full_first_author(self):
        first_author = self.authors.first()
        if first_author:
            return '{0} et al.'.format(first_author.print_full)
        else:
            return 'Unknown authors'

    @property
    def print_compact_first_author(self):
        first_author = self.authors.first()
        if first_author:
            return '{0} et al.'.format(first_author.print_compact)
        else:
            return 'Unknown authors'

    @property
    def print_json_ld_authors(self):
        authors = ['{0}, {1}'.format(a.last_name, a.first_name)
                   for a in self.authors.all()]
        if len(authors) > 1:
            return authors
        else:
            return authors[0]

    @property
    def print_date(self):
        if self.date_pp:
            return self.date_pp.strftime('%e %b %Y')
        elif self.date_ep:
            return '{date} (epub)'.format(date=self.date_ep.strftime('%e %b %Y'))
        elif self.date_fs:
            return '{date} (first seen)'.format(date=self.date_fs.strftime('%e %b %Y'))
        else:
            return 'Unknown Date'

    @property
    def print_publisher_name(self):
        if self.journal and self.journal.publisher:
            return self.journal.publisher.name
        else:
            return ''

    @property
    def print_month_year(self):
        if self.date_pp:
            return self.date_pp.strftime('%B %Y')
        elif self.date_ep:
            return '{date} (online)'.format(date=self.date_ep.strftime('%B %Y'))
        elif self.date_fs:
            return '{date} (first seen)'.format(date=self.date_fs.strftime('%B %Y'))
        else:
            return 'Unknown Date'

    @property
    def print_year(self):
        if self.date_pp:
            return self.date_pp.strftime('%Y')
        elif self.date_ep:
            return '{date} (online)'.format(date=self.date_ep.strftime('%Y'))
        elif self.date_fs:
            return '{date} (online)'.format(date=self.date_fs.strftime('%Y'))
        else:
            return 'Unknown Date'

    @property
    def print_first_seen(self):
        dates = [self.date_pp, self.date_ep, self.date_fs]
        dates = [d for d in dates if d]
        date = min(dates)
        return '{date}'.format(date=date.strftime('%e %b %Y'))

    @property
    def print_journal_title(self):
        if self.journal:
            return self.journal.title
        elif self.id_arx:
            return 'Arxiv'
        else:
            return 'Unknown journal'

    def print_journal_short_title(self):
        if self.journal:
            return self.journal.short_title or self.journal.title[:15] + '...'
        elif self.id_arx:
            return 'Arxiv'
        else:
            return 'Unknown journal'

    @property
    def print_ids(self):
        ids_str = ''
        # Loop through field starting with 'id_'
        for field in self._meta.fields:
            if field.name.startswith('id_') and getattr(self, field.name):
                if ids_str:
                    ids_str += ', '
                ids_str += '{id}: {value}'.format(
                    id=field.verbose_name,
                    value=getattr(self, field.name))
        return ids_str

    @property
    def print_clean_ids(self):
        """Not printing id_oth"""
        if self.id_doi:
            return 'DOI: {id}'.format(id=self.id_doi)
        elif self.id_arx:
            return 'Arxiv: {id}'.format(id=self.id_arx)
        elif self.id_pmi:
            return 'PMID: {id}'.format(id=self.id_pmi)
        elif self.id_pii:
            return 'PII: {id}'.format(id=self.id_pii)
        else:
            return ''

    @property
    def build_url(self):
        return 'https://%s%s' % (Site.objects.get_current(), self.get_absolute_url())

    def get_ids(self):
        """Return dictionary of paper ids"""
        ids = {}
        # Loop through field starting with 'id_'
        for field in self._meta.fields:
            if field.name.startswith('id_') and getattr(self, field.name):
                ids[field.name.split('id_')[1]] = getattr(self, field.name)
        return ids

    def count_ids(self):
        count = 0
        # Loop through field starting with 'id_'
        for field in self._meta.fields:
            if field.name.startswith('id_') and getattr(self, field.name):
                count += 1
        return count

    @staticmethod
    def detect_language(text):
        """ Detect language
        """
        lang_code = detect(text)
        return langcode_to_langpap(lang_code)

    def get_neighbors(self, time_span):
        from .tasks import get_neighbors_papers
        return get_neighbors_papers(self.id, time_span)

    def get_related_threads(self, time_span=-1):
        from .tasks import get_related_threads
        return get_related_threads(self.id, time_span)

    def state(self, user):
        if PaperUser.objects.filter(user=user, paper=self).exists():
            return PaperUser.objects.get(user=user, paper=self)
        else:
            return None

    def embed(self):
        from .tasks import embed_paper
        embed_paper(self.id)

    def get_doi_url(self):
        if self.id_doi:
            return 'http://doi.org/{doi}'.format(doi=self.id_doi)
        return None

    def __str__(self):
        return self.short_title


class AuthorPaper(TimeStampedModel):
    """Intermediate table for ranking authors in matches
    """
    author = models.ForeignKey(Author)
    paper = models.ForeignKey(Paper)
    position = models.IntegerField(default=1)

    def __str__(self):
        return '[{0}] {1}'.format(self.position, self.author.last_name)

    class Meta:
        ordering = ['position']


class CorpAuthorPaper(TimeStampedModel):

    corp_author = models.ForeignKey(CorpAuthor)
    paper = models.ForeignKey(Paper)

    def __str__(self):
        return '{0}'.format(self.corp_author.name)


class Stats(TimeStampedModel):

    nb_papers = models.IntegerField(default=0)

    nb_journals = models.IntegerField(default=0)

    nb_authors = models.IntegerField(default=0)

    nb_papers_last_week = models.IntegerField(default=0)

    nb_papers_last_two_weeks = models.IntegerField(default=0)

    nb_papers_last_month = models.IntegerField(default=0)

    nb_papers_last_year = models.IntegerField(default=0)

    def update(self):
        self.nb_papers = Paper.objects.count()
        self.nb_journals = Journal.objects.count()
        self.nb_authors = Author.objects.count()

        # Count matches by time range
        cursor = connection.cursor()
        query = "SELECT COUNT(*) FROM library_paper " \
                "WHERE date_co >= %s"

        d = timezone.now().date() - timezone.timedelta(days=7)
        cursor.execute(query, [d])
        self.nb_papers_last_week = cursor.fetchone()[0]

        d = timezone.now().date() - timezone.timedelta(days=14)
        cursor.execute(query, [d])
        self.nb_papers_last_two_weeks = cursor.fetchone()[0]

        d = timezone.now().date() - timezone.timedelta(days=30)
        cursor.execute(query, [d])
        self.nb_papers_last_month = cursor.fetchone()[0]

        d = timezone.now().date() - timezone.timedelta(days=365)
        cursor.execute(query, [d])
        self.nb_papers_last_year = cursor.fetchone()[0]

        self.save()

    def __str__(self):
        return self.created.strftime('%d %b %Y')

    def print(self):
        print('# matches: {matches}\n' \
               '# journals: {journals}\n' \
               '# authors: {authors}\n' \
               '# matches last week: {week}\n' \
               '# matches last 2 weeks: {two_weeks}\n' \
               '# matches last month: {month}\n' \
               '# matches last year: {year}\n'.format(
            papers=self.nb_papers,
            journals=self.nb_journals,
            authors=self.nb_authors,
            week=self.nb_papers_last_week,
            two_weeks=self.nb_papers_last_two_weeks,
            month=self.nb_papers_last_month,
            year=self.nb_papers_last_year,
        ))


class PaperUser(ModelDiffMixin, TimeStampedModel):
    # user
    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    # paper
    paper = models.ForeignKey(Paper)

    # Pinned or banned
    watch = models.PositiveIntegerField(null=True, default=None,
                                        choices=PAPER_WATCH,
                                        db_index=True)
    # Added or Trashed
    store = models.PositiveIntegerField(null=True, default=None,
                                        choices=PAPER_STORE,
                                        db_index=True)

    class Meta:
        unique_together = (('paper', 'user'),)

    def pin(self):
        self.watch = PAPER_PINNED
        self.save()

    def ban(self):
        self.watch = PAPER_BANNED
        self.save()

    def add(self, provider_id=None, info=None, orcid=False):
        """Add paper"""
        with transaction.atomic():
            if not provider_id:
                provider_id, info = self.user.lib.add_paper_on_provider(self.paper)
            self.user.lib.add_paper_on_etalia(self.paper, provider_id, info=info, orcid=orcid)
            self.store = PAPER_ADDED
            self.save()

    def trash(self):
        paper_provider_id = self.user.lib.userlib_paper\
            .get(paper=self.paper)\
            .paper_provider_id
        err = self.user.lib.trash_paper_on_provider(paper_provider_id)
        if not err:
            self.store = PAPER_TRASHED
            self.save()
            return None
        else:
            return err

    def save(self, **kwargs):
        if self.id:
            PaperUserHistory.objects.create(paperuser_id=self.id,
                                            difference=json.dumps(self.diff))
        super(PaperUser, self).save(**kwargs)


class PaperUserHistory(TimeStampedModel):

    paperuser = models.ForeignKey(PaperUser, related_name='history')

    difference = models.CharField(max_length=256, default='')

    date = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ('-created', )