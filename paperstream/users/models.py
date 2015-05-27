from django.db import models
from django.contrib.auth.models import AbstractBaseUser, \
    PermissionsMixin, BaseUserManager
from django.utils.translation import ugettext_lazy as _
from django.db.models import Q
from datetime import date
from django.core.mail import send_mail
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
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(**kwargs)
        user.is_superuser = True
        user.is_staff = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser, PermissionsMixin):

    # completely useless but required by python-social-auth
    username = models.CharField(_('username (UNUSED'), max_length=255,
                                blank=True, default='', db_index=True)

    email = models.EmailField(_('email address'), max_length=255,
                              unique=True, db_index=True)

    first_name = models.CharField(max_length=255, blank=True, default='')

    last_name = models.CharField(max_length=255, blank=True, default='')

    is_staff = models.BooleanField(_('staff status'), default=False,
        help_text=_('Designates whether the user can log into this admin '
                    'site.'))

    is_active = models.BooleanField(_('active'), default=True,
        help_text=_('Designates whether this user should be treated as '
                    'active. Unselect this instead of deleting accounts.'))

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


class UserLib(models.Model):
    """Library data of user"""

    user = models.OneToOneField(User, primary_key=True, related_name='userlib')

    affiliation = models.ForeignKey(Affiliation, null=True)

    papers = models.ManyToManyField(Paper, through='UserLibPaper')

    own_papers = models.ManyToManyField(Paper, related_name='own')

    journals = models.ManyToManyField(Journal, through='UserLibJournal')

    library_status = models.CharField(max_length=3, blank=True, default= '',
        choices=(('', 'Uninitialized'),
                 ('SYN', 'Synced'),
                 ('ING', 'Syncing')))

    feed_status = models.CharField(max_length=3, blank=True, default='',
        choices=(('', 'Uninitialized'),
                 ('SYN', 'Synced'),
                 ('ING', 'Syncing')))

    is_library_hooked = models.BooleanField(default=False)

    def count_papers(self):
        return self.papers.all().count()

    def count_journals(self):
        return self.papers.all().count()


class UserLibStats(models.Model):
    """Trace of user library/feed activity
    """
    userlib = models.ForeignKey(UserLib, related_name='stats')

    state = models.CharField(max_length=3,
                             choices=(('LOG', 'Log in'),
                                      ('LIB', 'Library syncing'),
                                      ('FEE', 'Feed syncing')))

    datetime = models.DateTimeField(null=False, auto_now_add=True)


class UserLibPaper(models.Model):
    """Intermediate model to user - paper
    """
    userlib = models.ForeignKey(UserLib)

    paper = models.ForeignKey(Paper)

    date_added = models.DateField(default=date(2000, 1, 1))

    class Meta:
        ordering = ['-date_added']

    def __str__(self):
        return '{0}@{1}'.format(self.paper.short_title(),
                                self.userlib.user.email)


class UserLibJournal(models.Model):
    """Intermediate model to user - journal
    """
    userlib = models.ForeignKey(UserLib)

    journal = models.ForeignKey(Journal)

    # number of paper link to that journal for this user
    papers_in_journal = models.IntegerField(default=0, null=False)

    # Score of the journal based on normalized distribution of paper per
    # journals
    score = models.FloatField(default=0.)

    def __str__(self):
        return '%s@%s' % (self.userlib.user.email, self.journal.short_title)

    def update_papers_in_journal(self):
        self.papers_in_journal = self.userlib.papers.filter(
            journal=self.journal).count()
        self.save()

    def update_score(self):
        self.update_papers_in_journal()
        total_papers_with_journal = \
            self.userlib.papers.filter(~Q(journal=None)).count()
        self.score = self.papers_in_journal / total_papers_with_journal
        self.save()



