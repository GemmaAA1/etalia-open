from nameparser import HumanName

from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from django.utils.translation import ugettext_lazy as _
from django.core.mail import send_mail
from django.conf import settings

from library.models import Paper, Journal
from nlp.models import Model
from core.constants import NLP_TIME_LAPSE_CHOICES
from feeds.constants import FEED_SCORING_CHOICES

from .validators import validate_first_name, validate_last_name
from core.models import TimeStampedModel


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


class UserManager(BaseUserManager):

    def create_user(self, email, password=None, **kwargs):
        email = self.normalize_email(email)
        user = self.model(email=email, is_active=True, **kwargs)
        user.set_password(password)
        user.save(using=self._db)
        UserLib.objects.create(user=user)
        UserStats.objects.log_user_init(user, '')
        UserSettings.objects.create(user=user)
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

    def set_state(self, state):
        if state in ['NON', 'IDL', 'ING']:
            self.state = state
            self.save()
            return self
        else:
            raise ValueError('Cannot set state. State value not allowed')


class UserLibPaper(TimeStampedModel):
    """Table - User/Paper relationship"""

    userlib = models.ForeignKey(UserLib)

    paper = models.ForeignKey(Paper)

    date_created = models.DateField(default=None, null=True)

    date_last_modified = models.DateField(default=None, null=True)

    authored = models.NullBooleanField(default=None, null=True, blank=True)

    starred = models.NullBooleanField(default=None, null=True, blank=True)

    scored = models.FloatField(default=0.)

    class Meta:
        ordering = ['-date_created']
        unique_together = [('userlib', 'paper')]

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title(),
                                self.userlib.user.email)


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

    def update_papers_in_journal(self):
        self.papers_in_journal = self.userlib.papers.filter(
            journal=self.journal).count()
        self.save()

    def __str__(self):
        return '%s@%s' % (self.userlib.user.email, self.journal.short_title)


class UserStatsManager(models.Manager):
    def log_lib_starts_sync(self, user, options=''):
        stats = self.model(user=user, state='LSS')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_lib_ends_sync(self, user, options=''):
        stats = self.model(user=user, state='LES')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_feed_starts_sync(self, user, options=''):
        stats = self.model(user=user, state='FSS')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_feed_ends_sync(self, user, options=''):
        stats = self.model(user=user, state='FES')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_init(self, user, options=''):
        stats = self.model(user=user, state='INI')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_email_valid(self, user, options=''):
        stats = self.model(user=user, state='EMA')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_log_in(self, user, options=''):
        stats = self.model(user=user, state='LIN')
        stats.options = options
        stats.save(using=self._db)
        return stats

    def log_user_log_out(self, user, options=''):
        stats = self.model(user=user, state='LOU')
        stats.options = options
        stats.save(using=self._db)
        return stats

class UserStats(models.Model):
    """Table - Log of UserLib and Feed activity"""

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


class UserSettingsManager(models.Manager):

    def create(self, **kwargs):
        # if nlp_model not defined, default to first nlp_model
        if 'model' not in kwargs:
            model = Model.objects.first()
            kwargs['model'] = model
        obj = self.model(**kwargs)
        obj.save(using=self._db)
        return obj


class UserSettings(TimeStampedModel):
    """Table - User Settings"""

    user = models.OneToOneField(settings.AUTH_USER_MODEL, primary_key=True,
                                related_name='settings')

    # NLP model to use
    model = models.ForeignKey(Model, verbose_name='NLP Model')

    # scoring method to use
    scoring_method = models.IntegerField(verbose_name='Scoring Algo',
                                         choices=FEED_SCORING_CHOICES,
                                         default=1)

    # in days
    time_lapse = models.IntegerField(default=61,
                                     choices=NLP_TIME_LAPSE_CHOICES,
                                     verbose_name='In the past for')

    objects = UserSettingsManager()

    def __str__(self):
        return self.user.email

