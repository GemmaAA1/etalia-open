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

from paperstream.library.models import Paper, Journal
from paperstream.nlp.models import Model
from paperstream.feeds.models import Stream
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


class UserLib(TimeStampedModel):
    """Table - User Library"""

    user = models.OneToOneField(User, primary_key=True, related_name='lib')

    papers = models.ManyToManyField(Paper, through='UserLibPaper')

    journals = models.ManyToManyField(Journal, through='UserLibJournal')

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

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title,
                                self.userlib.user.email)

    def print_created(self):
        return '{0:02d}-{1:02d}-{2:d}'.format(self.date_created.day,
                         self.date_created.month,
                         self.date_created.year)


class UserLibJournalManager(models.Manager):

    def add(self, **kwargs):
        obj, new = self.get_or_create(**kwargs)
        obj.papers_in_journal += 1
        obj.save()
        return obj


class UserLibJournal(TimeStampedModel):
    """Table - User/Journal relationship"""

    userlib = models.ForeignKey(UserLib)

    journal = models.ForeignKey(Journal)

    # number of paper link to that journal for this user
    papers_in_journal = models.IntegerField(default=0, null=False)

    objects = UserLibJournalManager()

    class Meta:
        unique_together = ('userlib', 'journal')
        ordering = ('-papers_in_journal', )

    def update_papers_in_journal(self):
        self.papers_in_journal = self.userlib.papers.filter(
            journal=self.journal).count()
        self.save()

    def __str__(self):
        return '%s@%s' % (self.userlib.user.email, self.journal.short_title)


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

    # in days
    stream_time_lapse = models.IntegerField(default=NLP_TIME_LAPSE_CHOICES[2][0],
                                            choices=NLP_TIME_LAPSE_CHOICES,
                                            verbose_name='Time range')

    # arbitrary units
    stream_narrowness = models.IntegerField(default=NLP_NARROWNESS_CHOICES[2][0],
                                            choices=NLP_NARROWNESS_CHOICES,
                                            verbose_name='Narrowness')

    # Trend settings
    # nlp model
    trend_model = models.ForeignKey(Model, verbose_name='NLP Model',
                                    related_name='trend_model')

    # scoring method to use
    trend_method = models.IntegerField(verbose_name='Method', default=0,
                                       choices=TREND_METHODS)

    # in days
    trend_time_lapse = models.IntegerField(default=NLP_TIME_LAPSE_CHOICES[2][0],
                                           choices=NLP_TIME_LAPSE_CHOICES,
                                           verbose_name='Time range')

    # arbitrary units
    trend_narrowness = models.IntegerField(default=NLP_NARROWNESS_CHOICES[2][0],
                                            choices=NLP_NARROWNESS_CHOICES,
                                            verbose_name='Narrowness')

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

    context_source = models.CharField(max_length=128)

    is_ticked = models.BooleanField(default=False)

    is_liked = models.BooleanField(default=False)

    class Meta:
        unique_together = [('user', 'paper'), ]

    def __str__(self):
        if self.is_liked:
            return '{user} liked {pk} from {context}'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.context_source)
        elif self.is_ticked:
            return '{user} disliked {pk} with {context}'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.context_source)
        else:
            return '{user}/{pk}/{context} has an issue'.format(
                user=self.user,
                pk=self.paper.id,
                context=self.context_source)


class FeedLayout(TimeStampedModel):
    """Store settings for stream display"""

    user = models.OneToOneField(User)

    stream_filter = JSONField(null=True)

    trend_filter = JSONField(null=True)

    library_filter = JSONField(null=True)



