from django.test import TestCase
from library.models import Paper, Journal, Author, AuthorPosition


class LibraryTestCase(TestCase):

    def journal_title(self, **kwargs):
        journal = Journal(title='Journal title', **kwargs)
        return journal

    def journal1(self):
        journal = Journal(id_issn='1053-8119', title='Journal title #1')
        return journal

    def journal2(self):
        journal = Journal(id_issn='0000-0019', title='Journal title #2')
        return journal

    def journal3(self):
        journal = Journal(id_issn='1476-4687', title='Journal title #3')
        return journal