from django.test import TestCase
from ..models import Paper, Journal, Author, Publisher, AuthorPaper
from django.core.exceptions import ValidationError
from django.utils import timezone

class PublisherModelTest(TestCase):

    def test_str_is_name(self):
        publisher = Publisher(name='Publisher #1')
        self.assertEqual(publisher.__str__(), 'Publisher #1')


class PaperModelTest(TestCase):

    def test_create_paper_empty(self):
        paper = Paper(title='Paper X')
        paper.full_clean()

    def test_get_absolute_url(self):
        paper = Paper(title='Paper X')
        paper.save()
        self.assertEqual(paper.get_absolute_url(),
                         '/library/paper/{0}/'.format(paper.pk, ))

    def test_paper_ordering(self):
        p1 = Paper.objects.create(title='Paper X', date_ep=timezone.now().date())
        p2 = Paper.objects.create(title='Paper Y',
            date_ep=(timezone.now() - timezone.timedelta(days=2)).date())
        self.assertEqual(list(Paper.objects.all()), [p1, p2])

    def test_paper_can_have_multiple_ids(self):
        Paper(id_doi='xxx', id_arx='yyy')

    def test_duplicate_doi_are_invalid(self):
        Paper.objects.create(id_doi='xxx')
        with self.assertRaises(ValidationError):
            Paper(id_doi='xxx').full_clean()

    def test_duplicate_doi_are_invalid_through_update(self):
        Paper.objects.create(id_doi='xxx')
        paper = Paper.objects.create(id_doi='yyy')
        paper.id_doi = 'xxx'
        with self.assertRaises(ValidationError):
            paper.full_clean()

    def test_duplicate_arx_are_invalid(self):
        Paper.objects.create(id_arx='xxx')
        with self.assertRaises(ValidationError):
            Paper(id_arx='xxx').full_clean()

    def test_duplicate_pmi_are_invalid(self):
        Paper.objects.create(id_pmi='xxx')
        with self.assertRaises(ValidationError):
            Paper(id_pmi='xxx').full_clean()

    def test_duplicate_oth_are_invalid(self):
        Paper.objects.create(id_oth='xxx')
        with self.assertRaises(ValidationError):
            Paper(id_oth='xxx').full_clean()

    def test_duplicate_none_identifiers_are_valid(self):
        Paper.objects.create()
        Paper.objects.create()

    def test_counts_id_defined(self):
        paper = Paper.objects.create(id_oth='xxx', id_arx='yyy')
        self.assertEqual(paper.count_ids(), 2)

    def test_print_ids(self):
        paper = Paper.objects.create(id_oth='xxx', id_arx='yyy',
                                     id_doi='0000-0019', id_pmi='pubmed id')
        disp_str = paper.print_ids
        self.assertIn('Arxiv: yyy', disp_str)
        self.assertIn('Other ID: xxx', disp_str)
        self.assertIn('DOI: 0000-0019', disp_str)
        self.assertIn('PMID: pubmed id', disp_str)


class JournalModelTest(TestCase):

    def setUp(self):
        Journal.objects.get_or_create(title='Journal X', id_issn='1053-8119')
        Journal.objects.get_or_create(title='Journal Z', id_issn='0000-0019')
        Journal.objects.get_or_create(title='Journal Y', id_issn='1476-4687')

    def test_get_absolute_url(self):
        journal = Journal.objects.first()
        self.assertEqual(journal.get_absolute_url(),
                         '/library/journal/{0}/'.format(journal.pk, ))

    def test_title_cannot_be_blank(self):
        with self.assertRaises(ValidationError):
            journal = Journal(title='')
            journal.full_clean()

    def test_duplicate_id_is_invalid(self):
        with self.assertRaises(ValidationError):
            journal = Journal.objects.first()
            journal_same = Journal(title='Journal Y', id_issn=journal.id_issn)
            journal_same.full_clean()

    def test_journals_ordered_by_title(self):
        titles = [journal.title for journal in Journal.objects.all()]
        self.assertEqual(titles, ['Journal X',
                                  'Journal Y',
                                  'Journal Z'])

class AuthorModelTest(TestCase):

    def test_empty_author_are_invalid(self):
        with self.assertRaises(ValidationError):
            author = Author(first_name='', last_name='')
            author.full_clean()

    def test_last_name_empty_is_invalid(self):
        with self.assertRaises(ValidationError):
            author = Author(last_name='', first_name='Bernard')
            author.full_clean()

    def test_default_email_is_empty(self):
        author = Author.objects.create(first_name='Bernard')
        self.assertEqual(author.email, '')

    def test_default_first_name_is_empty(self):
        author = Author.objects.create(last_name='Bernard')
        self.assertEqual(author.first_name, '')

    def test_str(self):
        author = Author(first_name='Bernard', last_name='Paul')
        self.assertEqual(author.__str__(), 'Bernard Paul')

    def test_compact_print(self):
        author = Author(first_name='Bernard', last_name='Le Normand')
        self.assertEqual(author.print_compact, 'Le Normand B.')
        author = Author(first_name='', last_name='Le Normand')
        self.assertEqual(author.print_compact, 'Le Normand')

    def test_print_full(self):
        author = Author(first_name='Patrick', last_name='Hernandez')
        self.assertEqual(author.print_full, 'Patrick Hernandez')
        author = Author(first_name='', last_name='Hernandez')
        self.assertEqual(author.print_full, 'Hernandez')


class PaperAuthorTest(TestCase):

    def setUp(self):
        Author.objects.create(first_name='Bernard', last_name='Weiner')
        Author.objects.create(first_name='', last_name='Christen')
        Author.objects.create(first_name='Louis Andre', last_name='Van Brook')
        Paper.objects.create(title='On the road again')

    def test_authors_ordering(self):
        authors = Author.objects.all()
        paper = Paper.objects.first()
        for i, author in enumerate(authors):
            AuthorPaper.objects.create(author=author,
                                          paper=paper,
                                          position=i+1)
        cs = [cp.author for cp in AuthorPaper.objects.filter(paper=paper)]
        self.assertEqual(cs, list(authors))

    #TODO: test display

class JournalPaperTest(TestCase):

    def test_paper_can_be_added_to_journal(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        Paper.objects.create(journal=journal)

    def test_counts_number_papers_in_journal(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        Paper.objects.create(journal=journal, id_doi='xxx')
        Paper.objects.create(journal=journal, id_doi='yyy')
        Paper.objects.create(journal=journal, id_doi='zzz')

        journal.count_papers()
        journal.save()
        journal = Journal.objects.get(id_issn='0000-0019')
        self.assertEqual(journal.paper_set.count(), 3)
        self.assertEqual(journal.lib_size, 3)

    def test_ordering_of_papers_in_journal_by_date(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        p1 = Paper.objects.create(id_doi='xxx', journal=journal,
                                  date_ep=timezone.datetime(2015, 1, 1))
        p2 = Paper.objects.create(id_doi='yyy', journal=journal,
                                  date_ep=timezone.datetime(2015, 1, 2))
        p3 = Paper.objects.create(id_doi='zzz', journal=journal,
                                  date_ep=timezone.datetime(2015, 1, 3))
        self.assertEqual(list(journal.paper_set.all()), [p3, p2, p1])