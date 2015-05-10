from django.db import models
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from stdnum import issn as issn_checker
from django.db.models import Q
from jsonfield import JSONField
import collections


class Publisher(models.Model):
    """Publisher group
    """
    # Group name
    name = models.CharField(max_length=200, unique=True)
    # base url
    url = models.URLField()

    def __str__(self):
        return self.name


class Journal(models.Model):
    """Periodicals
    """
    # Identifiers
    # TODO: define custom field for this ids
    id_issn = models.CharField(max_length=9, null=False, blank=True, default='',
                               db_index=True)
    id_eissn = models.CharField(max_length=9, null=False, blank=True, default='',
                                db_index=True)
    id_arx = models.CharField(max_length=32, null=False, blank=True, default='',
                                db_index=True)
    id_oth = models.CharField(max_length=32, null=False, blank=True, default='',
                                db_index=True)

    # periodical title
    title = models.CharField(max_length=200, blank=True, default='')
    # short title
    short_title = models.CharField(max_length=100, blank=True, default='')

    # publisher group
    publisher = models.ForeignKey(Publisher, null=True, default='', blank=True)
    # url
    url = models.URLField(blank=True, null=True, default='')
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

    def validate_unique_ids(self):
        if self.id_issn and \
                self._default_manager.filter(Q(id_issn=self.id_issn) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_issn already exists')
        if self.id_arx and \
                self._default_manager.filter(Q(id_arx=self.id_arx) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_arx already exists')
        if self.id_oth and \
                self._default_manager.filter(Q(id_oth=self.id_oth) &
                        ~Q(pk=self.pk)):
            raise ValidationError('id_oth already exists')

    def validate_unique(self, exclude=None):
        super(Journal, self).validate_unique(exclude=exclude)
        self.validate_unique_ids()

    # to force django to store None for empty string
    def clean_id_issn(self):
        if self.id_issn:
            issn_checker.validate(self.id_issn)

    def clean(self, *args, **kwargs):
        self.clean_id_issn()
        super(Journal, self).clean()

    def save(self, *args, **kwargs):
        self.validate_unique()
        self.clean_id_issn()
        super(Journal, self).save(*args, **kwargs)

    def __str__(self):
        if self.short_title:
            return self.short_title
        else:
            if len(self.title) > 30:
                return self.title[:30]+'...'
            else:
                return self.title[:30]

    def get_absolute_url(self):
        return reverse('view_journal', args=[self.id])

    def counts_paper(self):
        self.lib_size = len(self.paper_set.all())
        return self.lib_size


class Author(models.Model):
    """Creators (authors) of papers
    """
    # TODO: add affiliation of authors ?

    # first name
    first_name = models.CharField(max_length=100, blank=True, default='')
    # last name
    last_name = models.CharField(max_length=100, blank=True, default='')
    # email
    email = models.EmailField(max_length=254, blank=True, null=True, default='')

    def __str__(self):
        return self.first_name + ' ' + self.last_name

    class Meta:
        unique_together = ('first_name', 'last_name')

    def clean(self):
        if not self.first_name and not self.last_name:
            raise ValidationError('Either first_name or last_name must '
                                  'be defined')


class Paper(models.Model):
    """Scientific papers
    """

    # identifiers (uniqueness defined thereafter)
    id_doi = models.CharField(max_length=32, null=False, blank=True, default='',
                              db_index=True)
    id_arx = models.CharField(max_length=32, null=False, blank=True, default='',
                              db_index=True)
    id_pmi = models.CharField(max_length=32, null=False, blank=True, default='',
                              db_index=True)
    id_oth = models.CharField(max_length=32, null=False, blank=True, default='',
                              db_index=True)

    # article title
    title = models.CharField(max_length=500, blank=True, default='')
    # authors
    authors = models.ManyToManyField(Author, through='AuthorPosition')
    # abstract
    abstract = models.TextField(blank=True, default='')
    # journal
    journal = models.ForeignKey(Journal, null=True, blank=True, db_index=True)
    # volume
    volume = models.CharField(max_length=200, blank=True, default='')
    # issue
    issue = models.CharField(max_length=200, blank=True, default='')
    # page
    page = models.CharField(max_length=200, blank=True, default='')
    # date published
    date = models.DateField(null=True, blank=True)
    # url where seen
    url = models.URLField(blank=True, null=True, default='')
    # date added
    date_added = models.DateTimeField(auto_now_add=True)

    # Booleans
    # article in press
    is_aip = models.BooleanField(default=False)
    # pre-print
    is_pre_print = models.BooleanField(default=False)
    # valid if all field have been verified or paper from 'official source'
    is_valid = models.BooleanField(default=False)

    class Meta:
        ordering = ['-date']

    def validate_unique_ids(self):
        if self.id_doi and \
                self._default_manager.filter(Q(id_doi=self.id_doi) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_doi already exists')
        if self.id_arx and \
                self._default_manager.filter(Q(id_arx=self.id_arx) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_arx already exists')
        if self.id_pmi and \
                self._default_manager.filter(Q(id_pmi=self.id_pmi) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_pmi already exists')
        if self.id_oth and \
                self._default_manager.filter(Q(id_oth=self.id_oth) &
                        ~Q(pk=self.pk)).exists():
            raise ValidationError('id_oth already exists')

    def validate_unique(self, exclude=None):
        self.validate_unique_ids()
        super(Paper, self).validate_unique(exclude=exclude)

    def save(self, *args, **kwargs):
        self.validate_unique()
        super(Paper, self).save(*args, **kwargs)

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
        return reverse('view_paper', args=[self.id])


class AuthorPosition(models.Model):
    """Intermediate table for ranking authors in papers
    """
    author = models.ForeignKey(Author, null=False)
    paper = models.ForeignKey(Paper, null=False)
    position = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return '[{0}] {1}'.format(self.position, self.author.last_name)

    class Meta:
        ordering = ['position']


class NewPaper(models.Model):
    """New papers
    """
    paper = models.OneToOneField(Paper)

    lib_size = models.IntegerField(default=0)

    def __str__(self):
        return self.paper.short_title

    def counts_papers(self):
        self.lib_size = self.paper.all()
        self.save()
        return self.lib_size