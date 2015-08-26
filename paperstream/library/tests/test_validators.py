from django.test import TestCase
from django.core.exceptions import ValidationError

from ..validators import validate_issn, validate_author_names


class TestValidateISSN(TestCase):

    def test_issn_invalid_checksum(self):
        with self.assertRaises(ValidationError):
            validate_issn('0000-0001')

    def test_issn_invalid_length(self):
        with self.assertRaises(ValidationError):
            validate_issn('0000-001')

    def test_issn_invalid_format(self):
        with self.assertRaises(ValidationError):
            validate_issn('0000-000X')