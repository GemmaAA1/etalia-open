from django.test import TestCase
from library.models import Paper, Journal, Author, AuthorPosition
from django.core.exceptions import ValidationError
from django.utils import timezone
from stdnum.exceptions import InvalidChecksum, InvalidFormat, InvalidLength


class PaperModelTest(TestCase):

    def test_create_paper_empty(self):
        Paper.objects.create()

    def test_get_absolute_url(self):
        paper = Paper.objects.create()
        self.assertEqual(paper.get_absolute_url(),
                         '/library/paper/{0}/'.format(paper.pk, ))

    def test_paper_ordering(self):
        p1 = Paper.objects.create(date=timezone.now().date())
        p2 = Paper.objects.create(date=(timezone.now() -
                                        timezone.timedelta(days=2)).date())
        self.assertEqual(list(Paper.objects.all()), [p1, p2])

    def test_paper_can_have_multiple_ids(self):
        Paper.objects.create(id_doi='xxx', id_arx='yyy')

    def test_duplicate_doi_are_invalid(self):
        Paper.objects.create(id_doi='xxx')
        with self.assertRaises(ValidationError):
            Paper.objects.create(id_doi='xxx')

    def test_duplicate_arx_are_invalid(self):
        Paper.objects.create(id_arx='xxx')
        with self.assertRaises(ValidationError):
            Paper.objects.create(id_arx='xxx')

    def test_duplicate_pmi_are_invalid(self):
        Paper.objects.create(id_pmi='xxx')
        with self.assertRaises(ValidationError):
            Paper.objects.create(id_pmi='xxx')

    def test_duplicate_oth_are_invalid(self):
        Paper.objects.create(id_oth='xxx')
        with self.assertRaises(ValidationError):
            Paper.objects.create(id_oth='xxx')

    def test_duplicate_empty_identifiers_are_valid(self):
        Paper.objects.create(id_doi='')
        Paper.objects.create(id_doi='')

    def test_prints(self):
        paper = Paper.objects.create(title='On the road again')
        # should not raise
        paper.__str__()


class JournalModelTest(TestCase):

    def test_can_save_journal_with_valid_issn(self):
        Journal.objects.create(id_issn='1053-8119')

    def test_canNOT_save_journal_with_invalid_issn(self):
        journal = Journal(id_issn='0000-0001')
        with self.assertRaises(InvalidChecksum):
            journal.full_clean()
        journal = Journal(id_issn='0000-001')
        with self.assertRaises(InvalidLength):
            journal.full_clean()
        journal = Journal(id_issn='X000-0001')
        with self.assertRaises(InvalidFormat):
            journal.full_clean()

    def test_get_absolute_url(self):
        journal = Journal.objects.create(id_issn='1053-8119')
        self.assertEqual(journal.get_absolute_url(),
                         '/library/journal/{0}/'.format(journal.id, ))

    def test_duplicate_issn_are_invalid(self):
        Journal.objects.create(id_issn='1053-8119')
        with self.assertRaises(ValidationError):
            journal = Journal(id_issn='1053-8119')
            journal.full_clean()

    def test_journals_ordered_by_title(self):
        j1 = Journal.objects.create(id_issn='1476-4687',
                                    title='A good journal')
        j2 = Journal.objects.create(id_issn='0000-0019',
                                    title='Curated journal')
        j3 = Journal.objects.create(id_issn='1053-8119',
                                    title='Bad journal')
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
        p1 = Paper.objects.create()
        AuthorPosition.objects.create(author=c1, paper=p1, position=1)
        AuthorPosition.objects.create(author=c2, paper=p1, position=2)
        AuthorPosition.objects.create(author=c3, paper=p1, position=2)

        cs = [cp.author for cp in AuthorPosition.objects.filter(paper=p1)]

        self.assertEqual(cs, [c1, c2, c3])


class JournalPaperTest(TestCase):

    def test_paper_can_be_added_to_journal(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        Paper.objects.create(journal=journal)

    def test_counts_number_papers_in_journal(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        Paper.objects.create(journal=journal, id_doi='xxx')
        Paper.objects.create(journal=journal, id_doi='yyy')
        Paper.objects.create(journal=journal, id_doi='zzz')

        journal.counts_paper()
        journal = Journal.objects.get(id_issn='0000-0019')
        self.assertEqual(journal.paper_set.count(), 3)
        self.assertEqual(journal.lib_size, 3)

    def test_ordering_of_papers_in_journal_by_date(self):
        journal = Journal.objects.create(id_issn='0000-0019')
        p1 = Paper.objects.create(id_doi='xxx', journal=journal,
                                  date=timezone.datetime(2015, 1, 1))
        p2 = Paper.objects.create(id_doi='yyy', journal=journal,
                                  date=timezone.datetime(2015, 1, 2))
        p3 = Paper.objects.create(id_doi='zzz', journal=journal,
                                  date=timezone.datetime(2015, 1, 3))
        self.assertEqual(list(journal.paper_set.all()), [p3, p2, p1])