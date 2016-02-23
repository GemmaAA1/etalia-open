# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import collections

from nameparser import HumanName
from jsonfield import JSONField
from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from django.utils.translation import ugettext_lazy as _
from django.core.mail import send_mail
from django.conf import settings

from paperstream.library.models import Paper, Journal, Author
from paperstream.nlp.models import Model
from paperstream.feeds.models import Stream, Trend
from paperstream.feeds.constants import STREAM_METHODS, TREND_METHODS
from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES, \
    NLP_NARROWNESS_CHOICES, EMAIL_DIGEST_FREQUENCY_CHOICES
from paperstream.core.models import TimeStampedModel

from .validators import validate_first_name, validate_last_name
from .constants import INIT_STEPS


class Affiliation(TimeStampedModel):
    """Table for Affiliations"""
    # TODO: Implement prepopulate table

    department = models.CharField(max_length=200, blank=True, default='')

    institution = models.CharField(max_length=200, blank=True, default='')

    city = models.CharField(max_length=50, blank=True, default='')

    state = models.CharField(max_length=20, blank=True, default='')

    country = models.CharField(max_length=50, blank=True, default='')

    def __str__(self):
        print_fields = ['department', 'institution', 'city', 'state', 'country']
        print_aff_l = []
        for f in print_fields:
            if getattr(self, f):
                print_aff_l.append(getattr(self, f))
        return ', '.join(print_aff_l)

    class Meta:
        unique_together = ('department', 'institution', 'city', 'state',
                           'country')

    @property
    def print_affiliation(self):
        print_fields = ['department', 'institution', 'city', 'state', 'country']
        print_aff_l = []
        for f in print_fields:
            if getattr(self, f):
                print_aff_l.append(getattr(self, f))
        return ', '.join(print_aff_l)

    @property
    def is_empty(self):
        if not (self.department or
                    self.institution or
                    self.city or
                    self.state or
                    self.country):
            return True
        else:
            return False


class UserManager(BaseUserManager):

    def create_user(self, email, password=None, **kwargs):
        email = self.normalize_email(email)
        user = self.model(email=email, is_active=True, **kwargs)
        user.set_password(password)
        user.username = '{first} {last}'.format(first=user.first_name,
                                                last=user.last_name)
        user.save(using=self._db)
        # create user library
        UserLib.objects.create(user=user)
        # create user settings
        UserSettings.objects.create(user=user)
        # create user default (main) feed
        Stream.objects.create_main(user=user)
        # create user default (main) feed
        Trend.objects.create(user=user, name='main')
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(**kwargs)
        user.is_superuser = True
        user.is_staff = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser, PermissionsMixin):
    """Table - PaperStream User"""

    # completely useless but required by python-social-auth
    username = models.CharField(_('username (UNUSED)'), max_length=255,
                                blank=True, default='', db_index=True)

    email = models.EmailField(_('Email'), max_length=255,
                              unique=True, db_index=True)

    first_name = models.CharField(_('First Name'), max_length=255, blank=True, default='',
                                  validators=[validate_first_name])

    last_name = models.CharField(_('Last Name'), max_length=255, blank=True, default='',
                                 validators=[validate_last_name])

    title = models.CharField(max_length=32, blank=True, default='')

    position = models.CharField(max_length=64, blank=True, default='')

    is_staff = models.BooleanField(_('staff status'), default=False,
        help_text=_('Designates whether the user can log into this admin '
                    'site.'))

    is_active = models.BooleanField(_('active'), default=True,
        help_text=_('Designates whether this user should be treated as '
                    'active. Unselect this instead of deleting accounts.'))

    is_alpha = models.BooleanField(_('alpha'), default=True,
        help_text=_('Designates whether this user should be treated as '
                    'an early adopter user.'))

    init_step = models.CharField(max_length=3, default='NON', choices=INIT_STEPS,
        help_text=_('Tag where init user stands'))

    affiliation = models.ForeignKey(Affiliation, null=True, default=None)

    photo = models.ImageField(upload_to='photos', null=True)

    objects = UserManager()

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = ['']

    class Meta:
        ordering = ('email', )

    def __str__(self):
        return self.email

    def clean(self):
        self.clean_first_name()
        self.clean_last_name()

    def clean_first_name(self):
        name = self.first_name
        name = name.replace('.', ' ').strip()
        self.first_name = ' '.join([n.capitalize() for n in name.split(' ')])

    def clean_last_name(self):
        Hname = HumanName(' '.join([self.first_name, self.last_name]))
        Hname.capitalize()
        self.last_name = Hname.last

    @property
    def is_admin(self):
        return self.is_staff

    def get_short_name(self):
        if self.first_name:
            first_names_cap = ''.join([n[0] for n in self.first_name.split(' ')])
            return '{0} {1}'.format(self.last_name, first_names_cap)
        else:
            return '{0}'.format(self.last_name)

    def get_full_name(self):
        first_name = []
        for f in self.first_name.split(' '):
            if len(f) == 1:
                first_name.append('{0}.'.format(f))
            else:
                first_name.append(f)
        first = ' '.join(first_name)

        return '{0} {1}'.format(first, self.last_name)

    def get_absolute_url(self):
        return '/users/{0}'.format(self.pk)

    def email_user(self, subject, message, from_email=None, **kwargs):
        send_mail(subject, message, from_email, [self.email], **kwargs)

    def has_perm(self, perm, obj=None):
        return True

    def has_module_perms(self, app_label):
        return True

    def get_paper_state(self, paper_id):
        return dict(UserTaste.get_state(paper_id, self.id),
                    **UserLibPaper.get_state(paper_id, self.id))

    def get_counters(self):
        return dict(UserTaste.get_counters(self.id),
                    **UserLibPaper.get_counters(self.id))


class UserLib(TimeStampedModel):
    """Table - User Library"""

    user = models.OneToOneField(User, primary_key=True, related_name='lib')

    papers = models.ManyToManyField(Paper, through='UserLibPaper')

    journals = models.ManyToManyField(Journal, through='UserLibJournal')

    authors = models.ManyToManyField(Author, through='UserLibAuthor')

    state = models.CharField(max_length=3, blank=True, default='NON',
        choices=(('NON', 'Uninitialized'),
                 ('IDL', 'Idle'),
                 ('ING', 'Syncing')))

    @property
    def count_papers(self):
        return self.papers.all().count()

    @property
    def count_journals(self):
        return self.journals.all().count()

    @property
    def count_authors(self):
        return self.papers.values('authors').count()

    def set_state(self, state):
        if state in ['NON', 'IDL', 'ING']:
            self.state = state
            self.save()
            return self
        else:
            raise ValueError('Cannot set state. State value not allowed')

    def update_authors(self):
        authors = self.papers\
            .values_list('authors', flat=True)
        counter = collections.Counter(authors)
        for a, v in counter.items():
            ula, _ = UserLibAuthor.objects.get_or_create(userlib=self,
                                                         author_id=a)
            ula.occurrence = v
            ula.save()

    def get_authors_dist(self):
        return UserLibAuthor.objects\
                    .filter(userlib=self)\
                    .values_list('author', 'occurrence')

    def update_journals(self):
        journals = self.papers\
            .exclude(journal__isnull=True)\
            .values_list('journal', flat=True)
        counter = collections.Counter(journals)
        for j, v in counter.items():
            ula, _ = UserLibJournal.objects.get_or_create(userlib=self,
                                                         journal_id=j)
            ula.occurrence = v
            ula.save()

    def get_journals_dist(self):
        return UserLibJournal.objects\
                    .filter(userlib=self)\
                    .values_list('journal', 'occurrence')


class UserLibPaper(TimeStampedModel):
    """Table - User/Paper relationship"""

    userlib = models.ForeignKey(UserLib, related_name='userlib_paper')

    paper = models.ForeignKey(Paper, related_name='userlib_paper')

    date_created = models.DateField(default=None, null=True)

    date_last_modified = models.DateField(default=None, null=True)

    authored = models.NullBooleanField(default=None, null=True, blank=True)

    starred = models.NullBooleanField(default=None, null=True, blank=True)

    scored = models.FloatField(default=0.)

    paper_provider_id = models.CharField(max_length=64, default='')

    is_trashed = models.BooleanField(default=False)

    class Meta:
        ordering = ['-date_created']
        unique_together = [('userlib', 'paper')]

    def print_created(self):
        return '{0:02d}-{1:02d}-{2:d}'.format(self.date_created.day,
                         self.date_created.month,
                         self.date_created.year)

    @classmethod
    def get_state(cls, paper_id, user_id):
        try:
            ulp = cls.objects.get(paper_id=paper_id, userlib__user_id=user_id)
            return {'is_added': True, 'is_trashed': ulp.is_trashed}
        except cls.DoesNotExist:
            return {'is_added': False, 'is_trashed': False}

    @classmethod
    def get_counters(cls, user_id):
        try:
            objs = cls.objects.filter(userlib__user_id=user_id)\
                .values('is_trashed')
            return {
                'trash': len([obj for obj in objs if obj['is_trashed']]),
                'library': len([obj for obj in objs if not obj['is_trashed']]),
            }
        except cls.DoesNotExist:
            return {'trash': None, 'library': None}

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title,
                                self.userlib.user.email)


class UserLibJournalManager(models.Manager):

    def add(self, **kwargs):
        obj, new = self.get_or_create(**kwargs)
        obj.occurrence += 1
        obj.save()
        return obj


class UserLibJournal(TimeStampedModel):
    """Table - User/Journal relationship"""

    userlib = models.ForeignKey(UserLib)

    journal = models.ForeignKey(Journal)

    occurrence = models.IntegerField(default=0, null=False)

    objects = UserLibJournalManager()

    class Meta:
        unique_together = ('userlib', 'journal')
        ordering = ('-occurrence', )

    def update_occurrence(self):
        self.occurrence = self.userlib.papers.filter(
            journal=self.journal).count()
        self.save()

    def __str__(self):
        return '%s@%s' % (self.userlib.user.email, self.journal.short_title)


class UserLibAuthor(TimeStampedModel):

    userlib = models.ForeignKey(UserLib)

    author = models.ForeignKey(Author)

    occurrence = models.IntegerField(default=0, null=False)

    class Meta:
        unique_together = ('userlib', 'author')
        ordering = ('-occurrence', )

    def __str__(self):
        return '%s' % (self.author.last_name)

    def update_occurence(self):
        self.occurence = self.userlib.papers.filter(author=self.author).count()
        self.save()


class UserStatsManager(models.Manager):
    """"""
    def log_lib_starts_sync(self, user, options=''):
        stats = self.model(user=user, message='library start sync')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_lib_ends_sync(self, user, options=''):
        stats = self.model(user=user, message='library ends sync')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_stream_starts_sync(self, user, options=''):
        stats = self.model(user=user, message='stream start update')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_stream_ends_sync(self, user, options=''):
        stats = self.model(user=user, message='stream ends update')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_trend_starts_sync(self, user, options=''):
        stats = self.model(user=user, message='trend start update')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_trend_ends_sync(self, user, options=''):
        stats = self.model(user=user, message='trend ends update')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_log_in(self, user, options=''):
        stats = self.model(user=user, message='logging in')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_log_out(self, user, options=''):
        stats = self.model(user=user, message='logging out')
        stats.options = options
        stats.save(using=self._db)
        return stats


class UserStats(TimeStampedModel):
    """Table - Log of UserLib and Feed activity"""

    user = models.ForeignKey(User, related_name='stats')

    message = models.CharField(max_length=128)

    options = models.CharField(max_length=128, default='', blank=True)

    objects = UserStatsManager()


class UserSettingsManager(models.Manager):

    def create(self, **kwargs):
        # if nlp_model not defined, default to first nlp_model
        if 'model' not in kwargs:
            model = Model.objects.first()
            kwargs['stream_model'] = model
            kwargs['trend_model'] = model
        obj = self.model(**kwargs)
        obj.save(using=self._db)
        return obj


class UserSettings(TimeStampedModel):
    """Table - User Settings"""

    user = models.OneToOneField(settings.AUTH_USER_MODEL, primary_key=True,
                                related_name='settings')

    ##  Stream settings
    # NLP model to use
    stream_model = models.ForeignKey(Model, verbose_name='NLP Model',
                                     related_name='stream_model')

    # scoring method to use
    stream_method = models.IntegerField(verbose_name='Method', default=0,
                                        choices=STREAM_METHODS)

    # author weight
    stream_author_weight = models.FloatField(default=1.0,
                                             verbose_name='Author weight')

    # journal weight
    stream_journal_weight = models.FloatField(default=1.0,
                                              verbose_name='Journal weight')

    # vector weight
    stream_vector_weight = models.FloatField(default=1.0,
                                             verbose_name='Content weight')

    # DEPRECATED
    # in days
    stream_time_lapse = models.IntegerField(default=NLP_TIME_LAPSE_CHOICES[2][0],
                                            choices=NLP_TIME_LAPSE_CHOICES,
                                            verbose_name='Time range')
    # DEPRECATED
    # arbitrary units
    stream_narrowness = models.IntegerField(default=NLP_NARROWNESS_CHOICES[2][0],
                                            choices=NLP_NARROWNESS_CHOICES,
                                            verbose_name='Narrowness')

    # DEPRECATED
    # stream method arguments
    stream_method_args = JSONField(null=True, default=None, blank=True)

    # Trend settings
    # nlp model
    trend_model = models.ForeignKey(Model, verbose_name='NLP Model',
                                    related_name='trend_model')

    # scoring method to use
    trend_method = models.IntegerField(verbose_name='Method', default=0,
                                       choices=TREND_METHODS)
    # DEPRECATED
    # in days
    trend_time_lapse = models.IntegerField(default=NLP_TIME_LAPSE_CHOICES[2][0],
                                           choices=NLP_TIME_LAPSE_CHOICES,
                                           verbose_name='Time range')

    # DEPRECATED
    # arbitrary units
    trend_narrowness = models.IntegerField(default=NLP_NARROWNESS_CHOICES[2][0],
                                            choices=NLP_NARROWNESS_CHOICES,
                                            verbose_name='Narrowness')
    # DEPRECATED
    # stream method arguments
    trend_method_args = JSONField(null=True, default=None, blank=True)

    # doc vector weight
    trend_doc_weight = models.FloatField(default=1.0,
                                         verbose_name='Content weight')

    # altmetric vector weight
    trend_altmetric_weight = models.FloatField(default=1.0,
                                               verbose_name='Altmetric weight')

    # Email digest
    email_digest_frequency = models.IntegerField(
        default=EMAIL_DIGEST_FREQUENCY_CHOICES[0][0],
        choices=EMAIL_DIGEST_FREQUENCY_CHOICES,
        verbose_name='Email digest frequency')

    objects = UserSettingsManager()

    def __str__(self):
        return self.user.email


class UserTaste(TimeStampedModel):

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='tastes')

    paper = models.ForeignKey(Paper)

    source = models.CharField(max_length=128)

    is_banned = models.BooleanField(default=False)

    is_pinned = models.BooleanField(default=False)

    class Meta:
        unique_together = [('user', 'paper'), ]

    def __str__(self):
        if self.is_pinned:
            return '{user} pinned {pk} from {context}'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.source)
        elif self.is_banned:
            return '{user} banned {pk} with {context}'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.source)
        else:
            return '{user}/{pk}/{context} has an issue'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.source)

    @classmethod
    def get_state(cls, paper_id, user_id):
        try:
            ut = cls.objects.get(paper_id=paper_id, user_id=user_id)
            return {'is_pinned': ut.is_pinned, 'is_banned': ut.is_banned}
        except cls.DoesNotExist:
            return {'is_pinned': False, 'is_banned': False}

    @classmethod
    def get_counters(cls, user_id):
        try:
            objs = cls.objects.filter(user_id=user_id)\
                .values('is_pinned', 'is_banned')
            return {
                'pin': len([obj for obj in objs if obj['is_pinned']]),
                'ban': len([obj for obj in objs if obj['is_banned']])
            }
        except cls.DoesNotExist:
            return {'pin': None, 'ban': None}


class StreamLayout(TimeStampedModel):
    """Store settings for stream display"""

    user = models.OneToOneField(User)
    stream_filter = JSONField(null=True)


class TrendLayout(TimeStampedModel):
    """Store settings for trend display"""

    user = models.OneToOneField(User)
    trend_filter = JSONField(null=True)


class LibraryLayout(TimeStampedModel):
    """Store settings for library display"""

    user = models.OneToOneField(User)
    library_filter = JSONField(null=True)



