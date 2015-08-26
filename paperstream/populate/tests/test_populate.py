import os
import sys
from django.test import TestCase
from ..utils import populate_journal, populate_publisher
from paperstream.library.models import Journal, Publisher

JOU_LIST1 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'data_one_journal.csv')
JOU_LIST2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'data_two_journals.csv')
JOU_LIST3 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'data_arxiv_journals.csv')
PUB_LIST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'data_publisher.csv')


class PopulateBaseTest(TestCase):

    def pre_populate_publisher(self):
        populate_publisher(PUB_LIST, print_to=sys.stderr)


class PopulatePublisherTest(PopulateBaseTest):

    def test_publisher_can_be_save(self):
        records_added, errors = populate_publisher(PUB_LIST,
                                                   print_to=sys.stderr)
        self.assertEqual(records_added, Publisher.objects.count())
        self.assertEqual(len(errors), 0)

    def test_publisher_can_be_retrieved(self):
        populate_publisher(PUB_LIST, print_to=sys.stderr)
        pub_names = [pub.name for pub in Publisher.objects.all()]
        self.assertIn('Elsevier', list(pub_names))
        self.assertIn('Arxiv', list(pub_names))


class PopulateJournalTest(PopulateBaseTest):

    def test_through_errors_publisher_not_in_list(self):
        rec, errors = populate_journal(JOU_LIST3, print_to=sys.stderr)
        self.assertTrue(len(errors) > 0)

    def test_journal_can_be_save(self):
        # Pre populate publisher
        self.pre_populate_publisher()

        records_added, errors = populate_journal(JOU_LIST1, print_to=sys.stderr)
        self.assertEqual(records_added, Journal.objects.count())
        self.assertEqual(len(errors), 0)

    def test_journal_title_is_cleaned(self):
        # NB: LIST3 has a title in lower case
        # Pre populate publisher
        self.pre_populate_publisher()

        populate_journal(JOU_LIST3, print_to=sys.stderr)
        journals = Journal.objects.all()
        self.assertIn(journals[0].title, ['Statistics - Applications',
                                          'Statistics - Computation'])
        self.assertIn(journals[1].title, ['Statistics - Applications',
                                          'Statistics - Computation'])

    def test_journal_can_be_recalled(self):
        # Pre populate publisher
        self.pre_populate_publisher()

        populate_journal(JOU_LIST2, print_to=sys.stderr)
        jou_names = [jou.title for jou in Journal.objects.all()]
        self.assertIn('Neuroimage-Clinical', list(jou_names))
        self.assertIn('Neuroimage', list(jou_names))

    def test_journal_info_are_not_overwritten(self):
        # Pre populate publisher
        self.pre_populate_publisher()

        populate_journal(JOU_LIST2, print_to=sys.stderr)
        # Neuroimage publisher should not be deleted
        populate_journal(JOU_LIST1, print_to=sys.stderr)
        journal = Journal.objects.get(title='NeuroImage')
        publisher = Publisher.objects.get(name='Elsevier')
        self.assertEqual(journal.publisher, publisher)

class PopulateConsumerTest(PopulateBaseTest):

    def test_consumer_pubmed(self):
        pass

    def test_consumer_elsevier(self):
        pass

    def test_consumer_arxiv(self):
        pass