from unittest import skip

from django.test import TestCase
from django.contrib.auth import get_user_model
from django.core.exceptions import ValidationError

from library.models import Paper, Journal
from nlp.models import Model

from ..models import UserLib, UserStats, UserLibPaper, UserLibJournal, \
    UserSettings
from .base import LibDataTestCase

User = get_user_model()


class UserTest(TestCase):

    def test_user_can_be_instantiated(self):
        User(email='test@test.com')

    def test_user_must_have_email(self):
        with self.assertRaises(ValidationError):
            User().full_clean()

    def test_user_canNOT_have_same_email(self):
        user = User(email='test@test.com')
        user.save()
        user2 = User(email='test@test.com')
        with self.assertRaises(ValidationError):
            user2.full_clean()

    def test_email_must_be_valid(self):
        user = User(email='test.com')
        with self.assertRaises(ValidationError):
            user.full_clean()

    def test_default_user_is_not_staff_nor_admin(self):
        user = User(email='test@test.com')
        self.assertFalse(user.is_staff)
        self.assertFalse(user.is_admin)

# TODO: use python-nameparser instead
@skip
class UserTestNames(TestCase):

    def setUp(self):
        self.user = User(email='test@test.com',
                         first_name='James Douglas',
                         last_name='Morrison')
        self.user2 = User(email='test@test.com',
                          first_name='James d.',
                          last_name='Morrison')
        self.user3 = User(email='test@test.com',
                          first_name='j. D.',
                          last_name='Morrison')
        self.user4 = User(email='test@test.com',
                          first_name='J.D.',
                          last_name='Morrison')
        self.user5 = User(email='test@test.com',
                          first_name='JD',
                          last_name='Morrison-McCarthy')

    def test_user_can_display_short_name(self):
        self.assertEqual(self.user.get_short_name(),
                         'Morrison JD')
        self.assertEqual(self.user2.get_short_name(),
                         'Morrison JD')
        self.assertEqual(self.user3.get_short_name(),
                         'Morrison JD')

    def test_user_can_display_full_name(self):
        self.assertEqual(self.user.get_full_name(),
                         'James Douglas Morrison')
        self.assertEqual(self.user2.get_full_name(),
                         'James D. Morrison')
        self.assertEqual(self.user3.get_full_name(),
                         'J. D. Morrison')

class EmailUserTest(TestCase):
    #TODO: test email send address
    pass


class UserLibTest(LibDataTestCase):

    def setUp(self):
        super(UserLibTest, self).setUp()
        self.user = User(email='test@test.com')
        self.user.save()

    def test_user_lib_can_be_instantiated(self):
        ul = UserLib(user=self.user)
        ul.full_clean()

    def test_default_status_is_none(self):
        ul = UserLib(user=self.user)
        self.assertEqual(ul.state, 'NON')

    def test_lib_can_change_state(self):
        ul = UserLib.objects.create(user=self.user)
        ul.set_state('IDL')
        self.assertEqual(ul.state, 'IDL')

    def test_lib_can_add_journal(self):
        ul = UserLib.objects.create(user=self.user)
        ulj = UserLibJournal.objects.add(userlib=ul, journal=self.journal)
        self.assertTrue(ulj.papers_in_journal, 1)

    def test_lib_can_add_papers(self):
        ul = UserLib.objects.create(user=self.user)
        UserLibPaper.objects.create(userlib=ul, paper=self.paper)
        UserLibPaper.objects.create(userlib=ul, paper=self.paper2)
        self.assertTrue(ul.papers.count(), 2)


class UserStatsModelTest(TestCase):

    def setUp(self):
        super(UserStatsModelTest, self).setUp()
        self.user = User(email='test@test.com')
        self.user.save()

    def test_can_instantiate_user_stats(self):
        UserStats(user=self.user)

    def test_log_lib_starts_sync(self):
        UserStats.objects.log_lib_starts_sync(self.user)
        self.assertEqual(self.user.stats.first().state, 'LSS')

    def test_log_lib_ends_sync(self):
        UserStats.objects.log_lib_ends_sync(self.user, options='10')
        self.assertEqual(self.user.stats.first().state, 'LES')
        self.assertEqual(self.user.stats.first().options, '10')

    def test_log_feed_starts_sync(self):
        UserStats.objects.log_feed_starts_sync(self.user)
        self.assertEqual(self.user.stats.first().state, 'FSS')

    def test_log_feed_ends_sync(self):
        UserStats.objects.log_feed_ends_sync(self.user, options='10')
        self.assertEqual(self.user.stats.first().state, 'FES')
        self.assertEqual(self.user.stats.first().options, '10')

    def test_log_user_init(self):
        UserStats.objects.log_user_init(self.user)
        self.assertEqual(self.user.stats.first().state, 'INI')

    def test_log_user_email_valid(self):
        UserStats.objects.log_user_email_valid(self.user)
        self.assertEqual(self.user.stats.first().state, 'EMA')

    def log_user_log_in(self):
        UserStats.objects.log_user_log_in(self.user)
        self.assertEqual(self.user.stats.first().state, 'LIN')

    def log_user_log_out(self):
        UserStats.objects.log_user_log_out(self.user)
        self.assertEqual(self.user.stats.first().state, 'LOUT')


class UserSettingsTest(TestCase):

    def setUp(self):
        super(UserSettingsTest, self).setUp()
        self.user = User(email='test@test.com')
        self.user.save()
        self.model = Model.objects.create(name='test')

    def test_can_create_settings(self):
        UserSettings.objects.create(user=self.user, model=self.model)

    def test_can_create_settings_with_default_model(self):
        us = UserSettings.objects.create(user=self.user)
        self.assertEqual(us.model, self.model)

    def test_settings_time_lapse_is_in_choices(self):
        us = UserSettings.objects.create(user=self.user)
        us.time_lapse = 13214231412
        with self.assertRaises(ValidationError):
            us.full_clean()
