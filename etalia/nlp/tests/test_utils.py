# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from .base import NLPTestCase
from etalia.library.models import Paper, Journal
from ..utils import obj2tokens
from ..constants import FIELDS_FOR_MODEL

class TokenizePaperTest(NLPTestCase):

    def test_fields_must_be_kwargs(self):
        paper = Paper(title='Bla bla.', abstract='Hi. Hi, <p>hi</p> {mu}\n')
        with self.assertRaises(ValueError):
            obj2tokens(paper)

    def test_tokenize_paper(self):
        paper = Paper(title='Bla bla.', abstract='Hi. Hi, <p>hi</p> {mu}\n')
        default_fields = [field[0] for field in FIELDS_FOR_MODEL]
        tokens = obj2tokens(paper, fields=default_fields)
        self.assertEqual(tokens,
                         ['bla', 'bla', '.', 'hi', '.', 'hi', ',', 'hi',
                          '{', 'mu', '}', '.'])
