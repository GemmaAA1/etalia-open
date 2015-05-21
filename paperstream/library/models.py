from django.db import models
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from django.db.models import Q
from core.models import TimeStampedModel, NullableCharField
from .validators import validate_issn, validate_author_names
from .constants import LANGUAGES, PUBLISH_PERIODS, PAPER_TYPE, PUBLISH_STATUS
from .utils import langcode_to_langpap
from model_utils import Choices
from model_utils.fields import MonitorField, StatusField

from langdetect import detect

# Source from where Paper are created (tuple(set()) to make unique)

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

    # TODO: Test if db_index=True improve performance
    id_issn = NullableCharField(max_length=9, blank=True, null=True,
                                default=None, validators=[validate_issn],
                                unique=True, verbose_name='ISSN')
    id_eissn = NullableCharField(max_length=9, blank=True, null=True,
                                 default=None, validators=[validate_issn],
                                 unique=True, verbose_name='e-ISSN')
    id_arx = NullableCharField(max_length=32, blank=True, null=True,
                               default=None, unique=True,
                               verbose_name='Arxiv ID')
    id_oth = NullableCharField(max_length=32, blank=True, null=True,
                               default=None, unique=True,
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

    # Number of Paper in journal
    lib_size = models.IntegerField(default=0)

    # is flag
    is_trusted = models.BooleanField(default=False)

    class Meta:
        ordering = ['title']

    def __str__(self):
        if self.short_title:
            return self.short_title
        else:
            if len(self.title) > 30:
                return self.title[:30]+'...'
            else:
                return self.title[:30]

    def get_absolute_url(self):
        return reverse('library:journal', args=[self.id])

    def count_papers(self):
        self.lib_size = len(self.paper_set.all())
        return self.lib_size

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


class Author(TimeStampedModel):
    """Creators (authors) of papers
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
    """Scientific papers
    """

    # Published status
    publish_status = models.CharField(choices=PUBLISH_STATUS, blank=True,
                                      default='', max_length=20)
    publish_status_changed = MonitorField(monitor='publish_status')

    # Type of paper
    type = models.CharField(max_length=4, choices=PAPER_TYPE, blank=True,
                            default='')

    # TODO: Test if db_index=True improve performance
    # identifiers (uniqueness defined thereafter)
    id_doi = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='DOI')
    id_arx = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='Arxiv')
    id_pmi = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='PMID')
    # none unique because publisher dependent
    id_pii = NullableCharField(max_length=64, blank=True, default='',
                               null=True, verbose_name='PII')
    id_oth = NullableCharField(max_length=64, blank=True, default='',
                               null=True, unique=True, verbose_name='Other ID')
    # Title
    title = models.CharField(max_length=500, blank=False, default='')

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
    journal = models.ForeignKey(Journal, null=True, blank=True, default=None)
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
    date_ep = models.DateField(null=True, blank=True, default=None)
    # date of publication in issue journal
    date_pp = models.DateField(null=True, blank=True, default=None)
    # date of paper last revised (e.g. arxiv, or publisher with e.g grant#)
    date_lr = models.DateField(null=True, blank=True, default=None)

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
    is_trusted = models.BooleanField(default=False)

    def get_absolute_url(self):
        return reverse('library:paper', args=[self.pk])

    @property
    def short_title(self):
        return '{0}[...]'.format(self.title[:50])

    @property
    def print_compact_authors(self):
        authors = self.authors.all()
        authors_str = ''
        if authors:
            for author in authors:
                if authors_str:
                    authors_str += ', '
                authors_str += author.compact_disp
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
                authors_str += author.full_disp
            return authors_str
        else:
            return 'Unknown Authors'

    @property
    def print_full_first_author(self):
        first_author = self.authors.first()
        if first_author:
            return '{0} et al.'.format(first_author.full_disp)
        else:
            return 'Unknown authors'

    @property
    def print_compact_first_author(self):
        first_author = self.authors.first()
        if first_author:
            return '{0} et al.'.format(first_author.compact_disp)
        else:
            return 'Unknown authors'

    @property
    def print_date(self):
        if self.date:
            return self.date.strftime('%e %b %Y')
        else:
            return 'Unknown date'

    @property
    def print_journal_title(self):
        if self.journal:
            return self.journal.short_title or self.journal.title
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

    def __str__(self):
        return self.short_title


class AuthorPaper(TimeStampedModel):
    """Intermediate table for ranking authors in papers
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
