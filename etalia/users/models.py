# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from nameparser import HumanName
from mendeley.exception import MendeleyApiException

from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from django.utils.translation import ugettext_lazy as _
from django.core.mail import send_mail
from django.conf import settings
from django.db.models import F, Min
from django.utils import timezone


from etalia.library.models import Paper, Journal, Author, PaperUser
from etalia.feeds.constants import STREAM_METHODS, TREND_METHODS
from etalia.core.constants import EMAIL_DIGEST_FREQUENCY_CHOICES
from etalia.core.models import TimeStampedModel
from etalia.threads.models import Thread
from etalia.usersession.models import UserSession

from .validators import validate_first_name, validate_last_name
from .constants import INIT_STEPS, RELATIONSHIP_FOLLOWING, RELATIONSHIP_BLOCKED, \
    RELATIONSHIP_STATUSES, USERLIB_STATE_CHOICES, USERLIB_UNINITIALIZED, \
    USERLIB_NEED_REAUTH, USER_TYPES, USER_INDIVIDUAL


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
        # Stream.objects.create_main(user=user)
        # create user default (main) feed
        # Trend.objects.create(user=user, name='main')
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(**kwargs)
        user.is_superuser = True
        user.is_staff = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser, PermissionsMixin):
    """Table - etalia User"""

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

    init_step = models.CharField(max_length=3, default='NON', choices=INIT_STEPS,
        help_text=_('Tag where init user stands'))

    affiliation = models.ForeignKey(Affiliation, null=True, default=None)

    # photo = models.ImageField(upload_to='photos', null=True)
    #
    relationships = models.ManyToManyField('self', through='Relationship',
                                          symmetrical=False,
                                          related_name='related_to')

    # Define the type of user (individual, third-party, etc.)
    type = models.IntegerField(choices=USER_TYPES, default=USER_INDIVIDUAL,
                               null=False, blank=False, verbose_name='Type')

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

    # Relationships methods
    # ----------------------
    def add_relationship(self, user, status):
        relationship, created = Relationship.objects.get_or_create(
            from_user=self,
            to_user=user,
            status=status)
        return relationship, created

    def remove_relationship(self, user, status):
        Relationship.objects.filter(
            from_user=self,
            to_user=user,
            status=status).delete()
        return _, False

    def get_relationships(self, status):
        return self.relationships.filter(
            relation_to_users__status=status,
            relation_to_users__from_user=self)

    def get_related_to(self, status):
        return self.related_to.filter(
            relation_from_users__status=status,
            relation_from_users__to_user=self)

    def get_following(self):
        return self.get_relationships(RELATIONSHIP_FOLLOWING)

    def get_blocked(self):
        return self.get_relationships(RELATIONSHIP_BLOCKED)

    def get_followers(self):
        return self.get_related_to(RELATIONSHIP_FOLLOWING)

    def follow(self, user):
        return self.add_relationship(user, RELATIONSHIP_FOLLOWING)

    def block(self, user):
        return self.add_relationship(user, RELATIONSHIP_BLOCKED)

    def unfollow(self, user):
        return self.remove_relationship(user, RELATIONSHIP_FOLLOWING)

    def unblock(self, user):
        return self.remove_relationship(user, RELATIONSHIP_BLOCKED)

    @property
    def followers(self):
        return self.get_followers()

    @property
    def following(self):
        return self.get_following()

    @property
    def blocked(self):
        return self.get_blocked()


class UserLib(TimeStampedModel):
    """Table - User Library"""

    user = models.OneToOneField(User, primary_key=True, related_name='lib')

    papers = models.ManyToManyField(Paper, through='UserLibPaper')

    state = models.IntegerField(default=USERLIB_UNINITIALIZED,
                                choices=USERLIB_STATE_CHOICES)

    # date of first paper added in user library
    d_oldest = models.DateField(null=True, blank=True)

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
        valid_choices = [c[0] for c in USERLIB_STATE_CHOICES]
        if state in valid_choices:
            self.state = state
            self.save()
            return self
        else:
            raise ValueError('Cannot set state. State value not allowed')

    def get_d_oldest(self):
        # get date_created for papers
        created_date = self.userlib_paper.values('pk', 'date_created')
        # sort them
        created_date_s = sorted(created_date, key=lambda x: x['date_created'])
        return created_date_s[0]['date_created']

    def set_d_oldest(self):
        self.d_oldest = self.get_d_oldest()
        self.save()

    def get_session_backend(self):
        provider_name = self.user.social_auth.first().provider
        # get social
        social = self.user.social_auth.get(provider=provider_name)
        # get backend
        backend = social.get_backend_instance()
        # build session
        session = backend.get_session(social, self.user)
        return session, backend

    def add_paper_on_provider(self, paper):
        session, backend = self.get_session_backend()
        # add paper to provider
        err, id, info = backend.add_paper(session, paper)
        return id, info

    def add_paper_on_etalia(self, paper, provider_id, info=None):
        ulp, new = UserLibPaper.objects.get_or_create(userlib=self,
                                                      paper=paper)
        if info:
            ulp.date_created = info.get('created', None)
            ulp.last_date_modified = info.get('last_modified', None)
            ulp.authored = info.get('authored', None)
            ulp.starred = info.get('starred', None)
            ulp.scored = info.get('scored', 0.)
        ulp.paper_provider_id = provider_id
        ulp.save()

    def trash_paper_on_provider(self, provider_id):
        session, backend = self.get_session_backend()
        # remove paper from provider
        err = backend.trash_paper(session, provider_id)
        return err

    def update(self, full=False):
        try:
            session, backend = self.get_session_backend()
            # update lib
            backend.update_lib(self.user, session, full=full)
        except MendeleyApiException:
            # Put user on renew auth state
            self.set_state(USERLIB_NEED_REAUTH)
            UserSession.delete_user_sessions(self.user_id)

    def clear(self):
        ulps = UserLibPaper.objects.filter(userlib=self)
        pus = PaperUser.objects.filter(user=self.user)
        ulps.delete()
        pus.update(store=None)

    def reset(self):
        self.clear()
        self.update(full=True)


class UserLibPaper(TimeStampedModel):
    """Table - User/Paper relationship"""

    userlib = models.ForeignKey(UserLib, related_name='userlib_paper')

    paper = models.ForeignKey(Paper, related_name='userlib_paper')

    date_created = models.DateField(default=None, null=True, db_index=True)

    date_last_modified = models.DateField(default=None, null=True)

    authored = models.NullBooleanField(default=None, null=True, blank=True)

    starred = models.NullBooleanField(default=None, null=True, blank=True)

    scored = models.FloatField(default=0.)

    paper_provider_id = models.CharField(max_length=64, default='')

    # is_trashed = models.BooleanField(default=False)

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

    @property
    def is_trashed(self):
        # Deprecated with API
        if self.paper.paperuser_set.filter(user_id=self.id).exists():
            if self.paper.paperuser_set.filter(user_id=self.id).store == 2:
                return True
        return False


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


class UserSettings(TimeStampedModel):
    """Table - User Settings"""

    user = models.OneToOneField(settings.AUTH_USER_MODEL, primary_key=True,
                                related_name='settings')

    # User fingerprint

    # delta-time in MONTHS to roll back stream
    fingerprint_roll_back_deltatime = models.IntegerField(
        default=None,
        verbose_name='Roll-back time (months)',
        null=True,
        blank=True)

    # Stream settings
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
                                             verbose_name='Title/Abstract weight')

    # vector weight
    stream_score_threshold = models.FloatField(default=0.2,
                                               verbose_name='Specificity')


    # Trend settings

    # scoring method to use
    trend_method = models.IntegerField(verbose_name='Method', default=0,
                                       choices=TREND_METHODS)

    # doc vector weight
    trend_doc_weight = models.FloatField(default=0.2,
                                         verbose_name='Title/Abstract weight')

    # altmetric vector weight
    trend_altmetric_weight = models.FloatField(default=1,
                                               verbose_name='Altmetric weight')

    # vector weight
    trend_score_threshold = models.FloatField(default=0.1,
                                              verbose_name='Specificity')

    # Email digest
    email_digest_frequency = models.IntegerField(
        default=EMAIL_DIGEST_FREQUENCY_CHOICES[0][0],
        choices=EMAIL_DIGEST_FREQUENCY_CHOICES,
        verbose_name='Email digest frequency')

    def init_fingerprint_roll_back_deltatime(self):
        if not self.fingerprint_roll_back_deltatime:
            first_created = self.user.lib.userlib_paper.all()\
                .aggregate(min_date=Min(F('date_created')))['min_date']
            delta_months = (timezone.now().date() - first_created).days // \
                           (365 / 12) + 1
            self.fingerprint_roll_back_deltatime = delta_months
            self.save()
        return self.fingerprint_roll_back_deltatime

    def __str__(self):
        return self.user.email


class Relationship(TimeStampedModel):

    from_user = models.ForeignKey(User, related_name='relation_from_users')

    to_user = models.ForeignKey(User, related_name='relation_to_users')

    status = models.IntegerField(choices=RELATIONSHIP_STATUSES,
                                 default=RELATIONSHIP_FOLLOWING)

    def follow(self):
        self.status = RELATIONSHIP_FOLLOWING
        self.save()

    def block(self):
        self.status = RELATIONSHIP_BLOCKED
        self.save()

    class Meta:
        unique_together = ('from_user', 'to_user')


class UserInvited(TimeStampedModel):

    from_user = models.ForeignKey(settings.AUTH_USER_MODEL,
                                  null=True, blank=True, default=None)

    to_email = models.EmailField()
