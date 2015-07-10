from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from django.utils.translation import ugettext_lazy as _
from django.utils import timezone
from django.db.models import Q
from datetime import date
from django.core.mail import send_mail
from django.conf import settings
from django.core.validators import MinValueValidator, MaxValueValidator

from model_utils import fields

from library.models import Paper, Journal
from nlp.models import Model
from feeds.models import UserFeed
from feeds.constants import FEED_TIME_CHOICES

from .validators import validate_first_name, validate_last_name
from core.models import TimeStampedModel


class Affiliation(TimeStampedModel):
    """Affiliations
    """
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


class UserManager(BaseUserManager):

    def create_user(self, email, password=None, **kwargs):
        email = self.normalize_email(email)
        user = self.model(email=email, is_active=True, **kwargs)
        user.set_password(password)
        user.save(using=self._db)
        UserLib.objects.create(user=user)
        UserStats.objects.create_user_init(user, '')
        UserSettings.objects.create(user=user)
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(**kwargs)
        user.is_superuser = True
        user.is_staff = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser, PermissionsMixin):

    # completely useless but required by python-social-auth
    username = models.CharField(_('username (UNUSED)'), max_length=255,
                                blank=True, default='', db_index=True)

    email = models.EmailField(_('Email'), max_length=255,
                              unique=True, db_index=True)

    first_name = models.CharField(_('First Name'), max_length=255, blank=True, default='',
                                  validators=[validate_first_name])

    last_name = models.CharField(_('Last Name'), max_length=255, blank=True, default='',
                                 validators=[validate_last_name])

    is_staff = models.BooleanField(_('staff status'), default=False,
        help_text=_('Designates whether the user can log into this admin '
                    'site.'))

    is_active = models.BooleanField(_('active'), default=True,
        help_text=_('Designates whether this user should be treated as '
                    'active. Unselect this instead of deleting accounts.'))

    affiliation = models.ForeignKey(Affiliation, null=True, default=None)

    objects = UserManager()

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = ['']

    class Meta:
        ordering = ('email', )

    def __str__(self):
        return self.email

    @property
    def is_admin(self):
        return self.is_staff

    def get_short_name(self):
        if self.first_name:
            return '{0}. {1}'.format(self.first_name[0].capitalize(),
                                     self.last_name.capitalize())
        else:
            return '{0}'.format(self.last_name.capitalize())

    def get_full_name(self):
        return '{0} {1}'.format(self.first_name.capitalize(),
                                self.last_name.capitalize())

    def get_absolute_url(self):
        return '/users/{0}'.format(self.pk)

    def email_user(self, subject, message, from_email=None, **kwargs):
        send_mail(subject, message, from_email, [self.email], **kwargs)

    def has_perm(self, perm, obj=None):
        return True

    def has_module_perms(self, app_label):
        return True


class UserLib(TimeStampedModel):
    """Library data of user"""

    user = models.OneToOneField(User, primary_key=True, related_name='lib')

    papers = models.ManyToManyField(Paper, through='UserLibPaper')

    journals = models.ManyToManyField(Journal, through='UserLibJournal')

    status = models.CharField(max_length=3, blank=True, default='',
        choices=(('', 'Uninitialized'),
                 ('IDL', 'Idle'),
                 ('ING', 'Syncing')))

    @property
    def count_papers(self):
        return self.papers.all().count()

    @property
    def count_journals(self):
        return self.journals.all().count()

    def set_lib_syncing(self):
        self.status = 'ING'
        self.save()

    def set_lib_idle(self):
        self.status = 'IDL'
        self.save()


class UserStatsManager(models.Manager):
    def create_lib_starts_sync(self, user, options=''):
        stats = self.model(user=user, state='LSS')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_lib_ends_sync(self, user, options=''):
        stats = self.model(user=user, state='LES')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_feed_starts_sync(self, user, options=''):
        stats = self.model(user=user, state='FSS')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_feed_ends_sync(self, user, options=''):
        stats = self.model(user=user, state='FSS')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_user_init(self, user, options=''):
        stats = self.model(user=user, state='CRE')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_user_email_valid(self, user, options=''):
        stats = self.model(user=user, state='EMA')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_user_log_in(self, user, options=''):
        stats = self.model(user=user, state='LIN')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def create_user_log_out(self, user, options=''):
        stats = self.model(user=user, state='LOU')
        stats.options = options
        stats.save(using=self._db)
        return stats

class UserStats(models.Model):
    """Trace of user library/feed activity
    """
    user = models.ForeignKey(User, related_name='stats')

    state = models.CharField(max_length=3,
                             choices=(('LIN', 'Log in'),
                                      ('LOU', 'Log out'),
                                      ('LSS', 'Library starts syncing'),
                                      ('LES', 'Library ends syncing'),
                                      ('FSS', 'Feed starts sync'),
                                      ('FES', 'Feed ends sync'),
                                      ('EMA', 'Email validated'),
                                      ('CRE', 'Create user')))

    options = models.CharField(max_length=64, default='', blank=True)

    datetime = models.DateTimeField(null=False, auto_now_add=True)

    objects = UserStatsManager()


class UserLibPaper(TimeStampedModel):
    """Intermediate model to user - paper
    """
    userlib = models.ForeignKey(UserLib)

    paper = models.ForeignKey(Paper)

    date_created = models.DateField(default=date(2000, 1, 1))

    date_last_modified = models.DateField(default=date(2000, 1, 1))

    authored = models.NullBooleanField(default=None, null=True, blank=True)

    starred = models.NullBooleanField(default=None, null=True, blank=True)

    scored = models.FloatField(default=0.)

    class Meta:
        ordering = ['-date_created']

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title(),
                                self.userlib.user.email)

class UserLibJournal(TimeStampedModel):
    """Intermediate model to user - journal
    """
    userlib = models.ForeignKey(UserLib)

    journal = models.ForeignKey(Journal)

    # number of paper link to that journal for this user
    papers_in_journal = models.IntegerField(default=0, null=False)

    def __str__(self):
        return '%s@%s' % (self.userlib.user.email, self.journal.short_title)

    def update_papers_in_journal(self):
        self.papers_in_journal = self.userlib.papers.filter(
            journal=self.journal).count()
        self.save()


class UserSettingsManager(models.Manager):

    def create(self, **kwargs):
        # if nlp_model not defined, defaut to first nlp_model
        if 'model' not in kwargs:
            nlp_model = Model.objects.first()
            kwargs['model'] = nlp_model
        model = self.model(**kwargs)
        model.save(using=self._db)
        return model


class UserSettings(TimeStampedModel):

    user = models.OneToOneField(settings.AUTH_USER_MODEL, primary_key=True,
                                related_name='settings')

    # NLP model to use
    model = models.ForeignKey(Model, verbose_name='NLP Model')

    # in days
    time_lapse = models.IntegerField(default=7,
                                     choices=FEED_TIME_CHOICES,
                                     verbose_name='Feed from the past')

    objects = UserSettingsManager()

    def __str__(self):
        return self.user.email

