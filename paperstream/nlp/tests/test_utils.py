from django.test import TestCase
from library.models import Paper, Journal
from ..utils import paper2tokens
from ..constants import FIELDS_FOR_MODEL
from ..models import Journal, Model

class TokenizePaperTest(TestCase):
    def setUp(self):
        pass

    def test_fields_must_be_kwargs(self):
        paper = Paper(title='Bla bla.', abstract='Hi. Hi, <p>hi</p> {mu}\n')
        with self.assertRaises(ValueError):
            paper2tokens(paper)

    def test_tokenize_paper(self):
        paper = Paper(title='Bla bla.', abstract='Hi. Hi, <p>hi</p> {mu}\n')
        default_fields = [field[0] for field in FIELDS_FOR_MODEL]
        tokens = paper2tokens(paper, fields=default_fields)
        self.assertEqual(tokens,
                         ['bla', 'bla', '.', 'hi', '.', 'hi', ',', 'hi',
                          '{', 'mu', '}', '.'])
