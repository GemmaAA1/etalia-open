from django.db import models
from django.core.urlresolvers import reverse
from django.core.exceptions import ValidationError
from stdnum import issn as issn_checker


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
    # periodical title
    title = models.CharField(max_length=200, blank=True, default='')
    # short title
    short_title = models.CharField(max_length=100, blank=True, default='')
    # IDs
    # issn (print)
    issn = models.CharField(max_length=10, blank=True, default='')
    # issn (electronic or other)
    e_issn = models.CharField(max_length=10, blank=True, default='')
    # other type of ID, e.g. arxiv subcategory, conference proceedings
    ext_id = models.CharField(max_length=30, blank=True, default='')

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
        unique_together = ('issn', 'ext_id')
        ordering = ['title']

    def __str__(self):
        if self.short_title:
            return self.short_title
        else:
            if len(self.title) > 30:
                return self.title[:30]+'...'
            else:
                return self.title[:30]

    def clean(self):
        if self.ext_id and self.issn:
            raise ValidationError('Journal cannot have issn and ext_id IDs')
        if not self.ext_id and not self.issn:
            raise ValidationError('Journal must have an ID')

        # issn must be valid
        if self.issn:
            issn_checker.validate(self.issn)

    def get_absolute_url(self):
        return reverse('view_journal', args=[self.id])

    def counts_paper(self):
        self.lib_size = len(self.paper_set.all())
        self.save()
        return self.lib_size


class Author(models.Model):
    """Creators (authors) of papers
    """
    # TODO: add affiliation of creators ?

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
    # IDs
    # DOI
    doi = models.CharField(max_length=200, blank=True, default='')
    # External id, e.g Arxiv ID, Proceedings, ...
    ext_id = models.CharField(max_length=200, blank=True, default='')

    # article title
    title = models.CharField(max_length=500, blank=True, default='')
    # authors
    authors = models.ManyToManyField(Author, through='AuthorPosition')
    # abstract
    abstract = models.TextField(blank=True, default='')
    # journal
    journal = models.ForeignKey(Journal, null=True)
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
        unique_together = ('doi', 'ext_id')
        ordering = ['-date']

    def __str__(self):
        return self.doi or self.arxiv_id

    def clean(self):
        if self.ext_id and self.doi:
            raise ValidationError('Paper cannot have two different IDs')
        if not self.ext_id and not self.doi:
            raise ValidationError('Paper must have an ID')

    @property
    def short_title(self):
        return '{0}[...]'.format(self.title[:25])

    @property
    def creators_disp(self):
        cs = self.creators.all()
        if cs:
            # Init first guy
            ini = ''.join([fn[0] for fn in cs[0].first_name.split(' ')])
            creators_display = ' '.join((cs[0].last_name, ini+'.'))
            for c in cs[1:]:
                ini = ''.join([fn[0] for fn in c.first_name.split(' ')])
                creators_display = \
                    ', '.join((creators_display,
                               ' '.join((c.last_name, ini+'.'))))
            return creators_display
        else:
            return 'Unknown authors'

    @property
    def date_disp(self):
        if self.date:
            return self.date.strftime('%e %b %Y')
        else:
            return 'Unknown date'

    @property
    def title_disp(self):
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