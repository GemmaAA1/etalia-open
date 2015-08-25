from django.test import TestCase
from collections import namedtuple
from unittest.mock import patch

from ..pipeline import require_primary

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


class RequirePrimaryTest(TestCase):

    def setUp(self):
        self.basic_info = {'email': 'test@test.com',
                           'first_name': 'Paul',
                           'last_name': 'Touet'}
        session = {'session': {'basic_info': self.basic_info}}
        self.strategy = Struct(**session)
        self.details = {}
        self.user = Struct(**self.basic_info)

