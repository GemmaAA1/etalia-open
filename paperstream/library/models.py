from django.db import models
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from django.db.models import Q
from base.models import TimeStampedModel
from .validators import validate_id_issn, validate_id_eissn
import collections


class NullableCharField(models.CharField):
    description = "CharField that stores NULL but returns ''"
    __metaclass__ = models.SubfieldBase

    def to_python(self, value):
        if isinstance(value, models.EmailField):
            return value
        return value or ''

    def get_prep_value(self, value):
        return value or None


class Publisher(TimeStampedModel):
    """Publisher group
    """
    # Group name
    name = models.CharField(max_length=200, unique=True)
    # base url
    url = models.URLField()

    def __str__(self):
        return self.name


class Journal(TimeStampedModel):
    """Periodicals
    """
    # Identifiers
    # TODO: define custom field for this ids
    # TODO: Test if db_index=True improve performance
    id_issn = NullableCharField(max_length=9, blank=True, null=True,
                                default=None, validators=[validate_id_issn],
                                unique=True)
    id_eissn = NullableCharField(max_length=9, blank=True, null=True,
                                 default=None, validators=[validate_id_eissn],
                                 unique=True)
    id_arx = NullableCharField(max_length=32, blank=True, null=True,
                               default=None, unique=True)
    id_oth = NullableCharField(max_length=32, blank=True, null=True,
                               default=None, unique=True)

    # periodical title
    title = models.CharField(max_length=200)
    # short title
    short_title = models.CharField(max_length=100, blank=True)

    # publisher group
    publisher = models.ForeignKey(Publisher, null=True, default='', blank=True)
    # url
    url = models.URLField(blank=True, default='')
    # Scope
    scope = models.TextField(blank=True, max_length=1000, default='')

    # language
    language = models.CharField(max_length=200, blank=True, default='')

    # period
    PUBLISH_PERIOD = (('ANN', 'Annual'),
                      ('SEM', 'Semi-annual'),
                      ('TRI', 'Tri-annual'),
                      ('QUA', 'Quarterly'),
                      ('MON', 'Monthly'),
                      ('BIM', 'Bi-monthly'),
                      ('IRR', 'Irregular'))
    period = models.CharField(max_length=200, choices=PUBLISH_PERIOD,
                              default='IRR')

    # Number of Paper in journal
    lib_size = models.IntegerField(default=0)

    # is flag
    is_valid = models.BooleanField(default=False)

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

    def counts_papers(self):
        self.lib_size = len(self.paper_set.all())
        return self.lib_size

    def counts_ids(self):
        count = 0
        if self.id_issn:
            count += 1
        if self.id_i:
            count += 1
        if self.id_issn:
            count += 1
        if self.id_issn:
            count += 1
        return count

    def ids_disp(self):
        ids_build = []
        if self.id_issn:
            ids_build += 'ISSN: {0}'.format(self.id_issn)
        if self.id_eissn:
            if ids_build:
                ids_build += ', '
            ids_build += 'e-ISSN: {0}'.format(self.id_eissn)
        if self.id_arx:
            if ids_build:
                ids_build += ', '
            ids_build += 'Arxiv: {0}'.format(self.id_arx)
        if self.id_oth:
            if ids_build:
                ids_build += ', '
            ids_build += 'Other ID: {0}'.format(self.id_oth)
        return ids_build


class Author(TimeStampedModel):
    """Creators (authors) of papers
    """
    # TODO: add affiliation of authors ?

    # last name
    last_name = models.CharField(max_length=100)
    # first name
    first_name = models.CharField(max_length=100, blank=True, default='')
    # email
    email = models.EmailField(max_length=254, blank=True, default='')

    def __str__(self):
        return self.first_name + ' ' + self.last_name

    class Meta:
        unique_together = ('first_name', 'last_name')



class Paper(TimeStampedModel):
    """Scientific papers
    """
    # TODO: Test if db_index=True improve performance

    # identifiers (uniqueness defined thereafter)
    id_doi = NullableCharField(max_length=32, blank=True, default='',
                               null=True, unique=True)
    id_arx = NullableCharField(max_length=32, blank=True, default='',
                               null=True, unique=True)
    id_pmi = NullableCharField(max_length=32, blank=True, default='',
                               null=True, unique=True)
    id_oth = NullableCharField(max_length=32, blank=True, default='',
                               null=True, unique=True)
    # article title
    title = models.CharField(max_length=500)
    # authors
    authors = models.ManyToManyField(Author, through='AuthorPosition')
    # abstract
    abstract = models.TextField(blank=True, default='')
    # journal
    journal = models.ForeignKey(Journal, null=True, blank=True)
    # volume
    volume = models.CharField(max_length=200, blank=True, default='')
    # issue
    issue = models.CharField(max_length=200, blank=True, default='')
    # page
    page = models.CharField(max_length=200, blank=True, default='')
    # date published
    date = models.DateField(null=True, blank=True, )
    # url where seen
    url = models.URLField(blank=True, default='')

    # Booleans
    # article in press
    is_aip = models.BooleanField(default=False)
    # pre-print
    is_pre_print = models.BooleanField(default=False)
    # valid if all field have been verified or paper from 'official source'
    is_valid = models.BooleanField(default=False)

    class Meta:
        ordering = ['-date']

    def __str__(self):
        return self.short_title

    @property
    def short_title(self):
        return '{0}[...]'.format(self.title[:50])

    @property
    def authors_disp(self):
        cs = self.authors.all()
        if cs:
            # Init first guy
            ini = ''.join([fn[0] for fn in cs[0].first_name.split(' ')])
            authors_display = ' '.join((cs[0].last_name, ini+'.'))
            for c in cs[1:]:
                ini = ''.join([fn[0] for fn in c.first_name.split(' ')])
                authors_display = \
                    ', '.join((authors_display,
                               ' '.join((c.last_name, ini+'.'))))
            return authors_display
        else:
            return 'Unknown authors'

    @property
    def first_author_disp(self):
        first_author = self.authors.first()
        if first_author:
            return '{0}. {1} et al.'.format(
                first_author.last_name.capitalize(),
                first_author.first_name[0].capitalize())
        else:
            return 'Unknown authors'

    @property
    def date_disp(self):
        if self.date:
            return self.date.strftime('%e %b %Y')
        else:
            return 'Unknown date'

    @property
    def journal_title_disp(self):
        if self.journal:
            return self.journal.short_title.capitalize() \
                or self.journal.title.capitalize()
        else:
            return 'Unknown journal'

    def get_absolute_url(self):
        return reverse('library:paper', args=[self.pk])


class AuthorPosition(TimeStampedModel):
    """Intermediate table for ranking authors in papers
    """
    author = models.ForeignKey(Author)
    paper = models.ForeignKey(Paper)
    position = models.IntegerField(default=1)

    def __str__(self):
        return '[{0}] {1}'.format(self.position, self.author.last_name)

    class Meta:
        ordering = ['position']
