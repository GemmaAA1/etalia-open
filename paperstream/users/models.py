from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from datetime import date
from allauth.account.models import EmailAddress
from django.conf import settings
from django.core.validators import MinValueValidator, MaxValueValidator
from library.models import Paper, Journal
from model_utils import fields


class Affiliation(models.Model):
    """Affiliations
    """
    department = models.CharField(max_length=200, blank=True, default='')
    institution = models.CharField(max_length=200, blank=True, default='')
    city = models.CharField(max_length=50, blank=True, default='')
    state = models.CharField(max_length=20, blank=True, default='')
    country = models.CharField(max_length=50, blank=True, default='')

    def __str__(self):
        return self.name


class UserManager(BaseUserManager):
    def create_user(self, email, password, **kwargs):
        """
        Creates and saves a User with the given email and password.
        """
        if not email:
            raise ValueError('Users must have an email address')

        user = self.model(
            email=self.normalize_email(email),
            is_active=True,
            **kwargs
        )

        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_superuser(self, email, password, **kwargs):
        """
        Creates and saves a superuser with the given email, and password.
        """
        user = self.model(
            email=email,
            password=password,
            is_admin=True,
            is_superuser=True,
            is_active=True,
            **kwargs
        )
        user.is_admin = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser, PermissionsMixin):
    """user of the app
    """
    first_name = models.CharField(max_length=200, blank=True, default='')

    last_name = models.CharField(max_length=200, blank=True, default='')

    email = models.EmailField(
        verbose_name='Email address',
        max_length=255,
        unique=True,
        )

    # affiliation
    affiliation = models.ForeignKey(Affiliation, null=True)

    objects = UserManager()

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = ['']

    # Relationship to paper in library
    papers = models.ManyToManyField(Paper, through='UserPaper')

    # Relation to paper authored
    own_papers = models.ManyToManyField(Paper, related_name='own')

    lib_size = models.IntegerField(default=0)

    # Relationship to periodical
    journals = models.ManyToManyField(Journal, through='UserJournal')

    # Booleans
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)

    library_status = models.CharField(
        max_length=3,
        choices=(('SYN', 'Synced'), (('ING'), 'Syncing')))
    feed_status = models.CharField(
        max_length=3,
        choices=(('SYN', 'Synced'), (('ING'), 'Syncing')))

    is_library_hooked = models.BooleanField(default=False)

    def __str__(self):
        return self.email

    def has_perm(self, perm, obj=None):
        return True

    def has_module_perms(self, app_label):
        return True

    def account_verified(self):
        if self.is_authenticated:
            result = EmailAddress.objects.filter(email=self.email)
            if len(result):
                return result[0].verified
        return False

    def counts_papers(self):
        self.lib_size = len(self.papers.all())
        self.save()
        return self.lib_size

    @property
    def is_staff(self):
        return self.is_admin

    @property
    def lib_last_sync_str(self):
        return str(self.lib_last_sync)

    def get_short_name(self):
        if self.first_name:
            return '{0}. {1}'.format(self.first_name[0].capitalize(),
                                    self.last_name.capitalize())
        else:
            return '{0}'.format(self.last_name.capitalize())

    def get_full_name(self):
        return '{0} {1}'.format(self.first_name.capitalize(),
                                self.last_name.capitalize())


class UserStats(models.Model):

    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='stats')

    doing = models.CharField(max_length=3,
                             choices=(('LOG', 'Log in'),
                                      ('LIB', 'Library syncing'),
                                      ('FEE', 'Feed syncing')))

    when = models.DateTimeField(null=False)


class UserPaper(models.Model):
    """Intermediate model to link user to papers
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    paper = models.ForeignKey(Paper)
    date_added = models.DateField(default=date(2000, 1, 1))

    class Meta:
        ordering = ['-date_added']

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title(), self.user.email)


class UserJournal(models.Model):
    """Intermediate model to link user to journal
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL)
    journal = models.ForeignKey(Journal)

    # number of paper link to that journal for this user
    count = models.FloatField(default=0., null=False)

    # Last paper to be add from that journal
    last_addition = models.DateField(null=False, default=date(1900, 1, 1))

    # Score of the journal (TODO: to be defined, based on paper liked from
    # the journal ?)
    score = models.FloatField(default=1.)

    def __str__(self):
        return '%s@%s' % (self.user.email, self.journal.short_title)



