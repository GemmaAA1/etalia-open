# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.test import TestCase
from django.core.exceptions import ValidationError
from ..validators import validate_first_name, validate_last_name


class ValidateUserFirstNameTest(TestCase):

    def test_validate_first_name_raises_with_non_letter(self):
        with self.assertRaises(ValidationError):
            validate_first_name('$fir!')

    def test_validate_first_name_does_NOT_raise_with_odd_letter(self):
        validate_first_name('çaùmé')

    def test_validate_first_name_does_NOT_raise_with_quote(self):
        validate_first_name("van'much")

    def test_validate_first_name_can_have_space(self):
        validate_first_name('Tom Fr')

    def test_validate_first_name_cannot_have_period(self):
        with self.assertRaises(ValidationError):
            validate_first_name('Tom A.')


class ValidateUserLastNameTest(TestCase):

    def test_validate_last_name_raises_with_non_letter(self):
        with self.assertRaises(ValidationError):
            validate_last_name('$fir!')

    def test_validate_last_name_does_NOT_raise_with_odd_letter(self):
        validate_last_name('çaùmé')

    def test_validate_last_name_does_NOT_raise_with_quote(self):
        validate_last_name("van'much")

    def test_validate_last_name_can_have_space(self):
        validate_last_name('Tom Fr')

    def test_validate_last_name_cannot_have_period(self):
        with self.assertRaises(ValidationError):
            validate_last_name('Tom A.')
