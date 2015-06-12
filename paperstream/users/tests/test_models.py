from django.test import TestCase
from django.contrib.auth import get_user_model
from django.core.exceptions import ValidationError
from django.utils import timezone
from ..models import UserLib, UserStats
from library.models import Paper, Journal

User = get_user_model()


class TestCaseUser(TestCase):
    def setUp(self):
        User.objects.create_user(email='test@test.com')


class UserModelTest(TestCaseUser):

    def test_user_cannot_have_same_email(self):
        with self.assertRaises(ValidationError):
            User(email='test@test.com').full_clean()

    def test_email_must_be_valid(self):
        with self.assertRaises(ValidationError):
            User(email='test.com').full_clean()

    def test_create_user_create_lib_and_stats(self):
        user = User.objects.first()
        self.assertTrue(user.stats.exists)
        self.assertIsInstance(user.lib, UserLib)

    def test_default_user_is_not_staff(self):
        user = User.objects.first()
        self.assertFalse(user.is_staff)


class UserStatsModelTest(TestCaseUser):

    def test_create_user_row(self):
        user = User.objects.first()
        self.assertEqual(user.stats.first().state, 'CRE')

    def test_create_lib_stat_row(self):
        user = User.objects.first()
        UserStats.objects.create_lib_row(user, 10)
        self.assertEqual(user.stats.last().number_papers, 10)
        self.assertEqual(user.stats.last().state, 'LIB')

    def test_create_feed_stat_row(self):
        user = User.objects.first()
        UserStats.objects.create_feed_row(user, 10)
        self.assertEqual(user.stats.last().number_papers, 10)
        self.assertEqual(user.stats.last().state, 'FEE')

    def test_create_log_in_stat_row(self):
        user = User.objects.first()
        UserStats.objects.create_user_log_in_row(user)
        self.assertEqual(user.stats.last().state, 'LIN')

    def test_create_log_out_stat_row(self):
        user = User.objects.first()
        UserStats.objects.create_user_log_out_row(user)
        self.assertEqual(user.stats.last().state, 'LOU')


class UserLibModelTest(TestCaseUser):

    def test_change_lib_status_to_idle(self):
        user = User.objects.first()
        user.lib.set_lib_idle()
        self.assertEqual(user.lib.status, 'IDL')

    def test_change_lib_status_to_syncing(self):
        user = User.objects.first()
        user.lib.set_lib_syncing()
        self.assertEqual(user.lib.status, 'ING')


