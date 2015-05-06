from django.test import TestCase
from library.models import Paper, Journal, Author, AuthorPosition
from django.core.exceptions import ValidationError
from django.utils import timezone
from stdnum.exceptions import InvalidChecksum, InvalidFormat, InvalidLength

# Create your tests here.


class PaperModelTest(TestCase):
    def test_get_absolute_url(self):
        paper = Paper.objects.create(doi='xxx')
        self.assertEqual(paper.get_absolute_url(),
                         '/library/papers/{0}'.format(paper.id, ))

    def test_paper_ordering(self):
        paper1 = Paper.objects.create(doi='xxx',
                                      date=timezone.now().date())
        paper2 = Paper.objects.create(doi='yyy',
                                      date=(timezone.now() -
                                            timezone.timedelta(days=2)).date())
        self.assertEqual(
            list(Paper.objects.all()),
            [paper1, paper2])

    def test_empty_doi_and_empty_ext_id_are_invalid(self):
        paper = Paper.objects.create(doi='')
        with self.assertRaises(ValidationError):
            paper.ext_id = ''
            paper.full_clean()

    def test_duplicate_doi_are_invalid(self):
        Paper.objects.create(doi='xxx')
        with self.assertRaises(ValidationError):
            paper = Paper(doi='xxx')
            paper.full_clean()

    def test_duplicate_ext_id_are_invalid(self):
        Paper.objects.create(ext_id='xxx')
        with self.assertRaises(ValidationError):
            paper = Paper(ext_id='xxx')
            paper.full_clean()

    def test_paper_cannot_have_multiple_id(self):
        paper = Paper.objects.create(doi='xxx')
        with self.assertRaises(ValidationError):
            paper.ext_id = 'yyy'
            paper.full_clean()

    def test_paper_default_title_is_empty(self):
        paper = Paper.objects.create(doi='xxx')
        self.assertEqual(paper.title, '')


class JournalModelTest(TestCase):
    def test_get_absolute_url(self):
        journal = Journal.objects.create(issn='1053-8119')
        self.assertEqual(journal.get_absolute_url(),
                         '/library/journals/{0}'.format(journal.id, ))

    def test_empty_issn_and_empty_ext_id_are_invalid(self):
        journal = Journal.objects.create(issn='')
        with self.assertRaises(ValidationError):
            journal.ext_id = ''
            journal.full_clean()

    def test_duplicate_issn_are_invalid(self):
        Journal.objects.create(issn='1053-8119')
        with self.assertRaises(ValidationError):
            journal = Journal(issn='1053-8119')
            journal.full_clean()

    def test_duplicate_ext_id_are_invalid(self):
        Journal.objects.create(ext_id='xxx')
        with self.assertRaises(ValidationError):
            journal = Journal(ext_id='xxx')
            journal.full_clean()

    def test_journal_cannot_have_multiple_id(self):
        journal = Journal.objects.create(issn='1053-8119')
        with self.assertRaises(ValidationError):
            journal.ext_id = 'yyy'
            journal.full_clean()

    def test_journal_issn_must_be_valid(self):
        journal = Journal(issn='0000-0001')
        with self.assertRaises(InvalidChecksum):
            journal.full_clean()
        journal = Journal.objects.create(issn='0000-001')
        with self.assertRaises(InvalidLength):
            journal.full_clean()
        journal = Journal.objects.create(issn='X000-0001')
        with self.assertRaises(InvalidFormat):
            journal.full_clean()

    def test_journals_ordered_by_title(self):
        j1 = Journal.objects.create(issn='0000-0001', title='A good journal')
        j2 = Journal.objects.create(issn='0000-0021', title='Curated journal')
        j3 = Journal.objects.create(issn='1053-8119', title='Bad journal')
        js = Journal.objects.all()
        self.assertEqual(list(js), [j1, j3, j2])

class AuthorModelTest(TestCase):

    def test_empty_author_are_invalid(self):
        with self.assertRaises(ValidationError):
            author = Author(first_name='', last_name='')
            author.full_clean()

    def test_default_email_is_empty(self):
        author = Author.objects.create(first_name='Bernard')
        self.assertEqual(author.email, '')

    def test_default_first_name_is_empty(self):
        author = Author.objects.create(last_name='Bernard')
        self.assertEqual(author.first_name, '')


class PaperAuthorTest(TestCase):

    def test_authors_ordering(self):
        c1 = Author.objects.create(first_name='Bernard')
        c2 = Author.objects.create(first_name='Paul')
        c3 = Author.objects.create(first_name='Louis')
        p1 = Paper.objects.create(doi='xxx')
        AuthorPosition.objects.create(author=c1, paper=p1, position=1)
        AuthorPosition.objects.create(author=c2, paper=p1, position=2)
        AuthorPosition.objects.create(author=c3, paper=p1, position=2)

        cs = [cp.author for cp in AuthorPosition.objects.filter(paper=p1)]

        self.assertEqual(cs, [c1, c2, c3])


class JournalPaperTest(TestCase):

    def test_paper_can_be_added_to_journal(self):
        journal = Journal.objects.create(issn='0000-0000')
        paper = Paper.objects.create(doi='xxx', journal=journal)

    def test_counts_number_papers_in_journal(self):
        journal = Journal.objects.create(issn='0000-0000')
        Paper.objects.create(doi='xxx', journal=journal)
        Paper.objects.create(doi='yyy', journal=journal)
        Paper.objects.create(doi='zzz', journal=journal)
        journal.counts_paper()

        journal = Journal.objects.get(issn='0000-0000')
        self.assertEqual(journal.paper_set.count(), 3)
        self.assertEqual(journal.lib_size, 3)

    def test_ordering_of_papers_in_journal(self):
        journal = Journal.objects.create(issn='0000-0000')
        p1 = Paper.objects.create(doi='xxx', journal=journal,
                                  date=timezone.datetime(2015, 1, 1))
        p2 = Paper.objects.create(doi='yyy', journal=journal,
                                  date=timezone.datetime(2015, 1, 2))
        p3 = Paper.objects.create(doi='zzz', journal=journal,
                                  date=timezone.datetime(2015, 1, 3))

        self.assertEqual(list(journal.paper_set.all()), [p3, p2, p1])