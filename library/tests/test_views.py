from django.test import TestCase
from library.models import Paper, Journal, Author, Publisher, AuthorPosition
from django.core.exceptions import ValidationError


class LibraryView(TestCase):

    def test_library_page_render_library_template(self):
        response = self.client.get('/library/')
        self.assertTemplateUsed(response, 'library.html')


class JournalsView(TestCase):
    pass


class JournalView(TestCase):
    pass