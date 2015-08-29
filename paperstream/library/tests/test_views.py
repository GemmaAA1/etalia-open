# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from django.test import TestCase
from django.utils import timezone

from ..models import Paper, Journal, Author, Publisher, AuthorPaper
from . import gen_issn


class LibraryViewTest(TestCase):

    def test_library_page_render_library_template(self):
        response = self.client.get('/library/')
        self.assertTemplateUsed(response, 'library/library.html')


class JournalsViewTest(TestCase):

    def test_journals_page_render_journals_template(self):
        response = self.client.get('/library/journals/')
        self.assertTemplateUsed(response, 'library/journals.html')

    def test_view_passes_journals_list_to_template(self):
        Journal.objects.create(id_issn=gen_issn(),
                               title='On the Stack')
        Journal.objects.create(id_issn=gen_issn(),
                               title='In the Stack')
        journal_list = list(Journal.objects.all())
        response = self.client.get('/library/journals/')
        self.assertEqual(list(response.context['journal_list']), journal_list)


class JournalViewTest(TestCase):

    def test_journal_page_render_journal_template(self):
        journal = Journal.objects.create(id_issn=gen_issn(),
                                         title='On the Stack')
        response = self.client.get('/library/journal/{0}/'.format(journal.pk))
        self.assertTemplateUsed(response, 'library/journal.html')

    def test_passes_journal_info_to_template(self):
        journal = Journal.objects.create(id_issn=gen_issn(),
                                         title='On the Stack')
        response = self.client.get('/library/journal/{0}/'.format(journal.pk))
        self.assertEqual(response.context['journal'], journal)

    def test_passes_journal_info_and_paper_list_to_template(self):
        journal = Journal.objects.create(id_issn=gen_issn(),
                                         title='On the Stack')
        paper1 = Paper.objects.create(id_doi='xxx', journal=journal,
                                      title='A short story',
                                      date_ep=timezone.datetime(2015, 1, 1).date())
        paper2 = Paper.objects.create(id_doi='yyy', journal=journal,
                                      title='A long story',
                                      date_ep=timezone.datetime(2014, 1, 1).date())
        author1 = Author.objects.create(first_name='John', last_name='Crane')
        author2 = Author.objects.create(first_name='Keith', last_name='Haring')
        AuthorPaper.objects.create(paper=paper1, author=author1, position=0)
        AuthorPaper.objects.create(paper=paper2, author=author2, position=0)

        response = self.client.get('/library/journal/{0}/'.format(journal.pk))

        paper_list = list(Paper.objects.all())

        self.assertIn(response.context['paper_list'][0], paper_list)
        self.assertIn(response.context['paper_list'][1], paper_list)
