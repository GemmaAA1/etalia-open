import json
import sys
if sys.version_info < (3, 3):
    from mock import patch
else:
    from unittest.mock import patch

from django.test import TestCase
from django.core.exceptions import ValidationError

from paperstream.library.models import Journal, Paper
from ..models import ConsumerPubmed, ConsumerJournal, ConsumerElsevier, \
    ConsumerArxiv


class ConsumerTest(TestCase):

    def setUp(self):
        journal = Journal.objects.create(title='Journal of the consumers',
                                         id_issn='0000-0019')
        consumer = ConsumerPubmed.objects.create()
        consumer.day0 = 20
        consumer.save()
        cj = consumer.add_journal(journal)

    @staticmethod
    def return_utils():
        consumer = ConsumerPubmed.objects.first()
        journal = Journal.objects.first()
        cj = ConsumerJournal.objects.first()
        return consumer, journal, cj

    def test_consumer_have_unique_name(self):
        ConsumerPubmed.objects.create(name='consumer #1')
        with self.assertRaises(ValidationError):
            ConsumerPubmed(name='consumer #1').full_clean()

    def test_can_activate_journal(self):
        consumer, journal, cj = self.return_utils()
        consumer.activate_journal(journal)
        cj = consumer.consumerjournal_set.get(journal=journal)
        self.assertEqual(cj.status, 'idle')

    def test_can_deactivate_journal(self):
        consumer, journal, cj = self.return_utils()
        consumer.activate_journal(journal)
        consumer.deactivate_journal(journal)
        self.assertEqual(cj.status, 'inactive')

    def test_cannot_activate_if_error_status(self):
        consumer, journal, cj = self.return_utils()
        cj.status = 'error'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.activate_journal(journal)

    def test_cannot_deactivate_if_doing_stuff(self):
        consumer, journal, cj = self.return_utils()
        cj.status = 'consuming'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.deactivate_journal(journal)
        cj.status = 'in_queue'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.deactivate_journal(journal)


class ConsumerPubmedTest(TestCase):

    def setUp(self):
        journal, _ = Journal.objects.get_or_create(id_issn='1053-8119')
        consumer, _ = ConsumerPubmed.objects.get_or_create(name='pub cons #1')
        consumer.day0 = 10
        consumer.add_journal(journal)
        consumer.activate_journal(journal)
        consumer.save()

        with open('paperstream/consumers/tests/pubmed_sample.json') as file:
            self.entries = json.load(file)

    def test_add_entry(self):
        consumer = ConsumerPubmed.objects.first()
        journal = Journal.objects.first()
        item = consumer.parser.parse(self.entries[0])
        consumer.add_entry(item, journal)
        self.assertEqual(Paper.objects.all().count(), 1)

    @patch('paperstream.consumers.models.ConsumerPubmed.consume_journal')
    def test_consume_update_journal(self, mock_consume_journal):
        # setup mock object
        mock_consume_journal.return_value = self.entries, True

        journal = Journal.objects.first()
        consumer = ConsumerPubmed.objects.first()
        consumer.populate_journal(journal.id)

        self.assertIsNotNone(consumer.consumerjournal_set.first().last_date_cons)
        self.assertTrue(
            consumer.consumerjournal_set.first().last_number_papers_recorded > 0)
        self.assertTrue(
            consumer.consumerjournal_set.first().last_number_papers_fetched > 0)

    @patch('paperstream.consumers.models.ConsumerPubmed.consume_journal')
    def test_populate_journal(self, mock_consume_journal):
        # setup mock object
        mock_consume_journal.return_value = self.entries, True

        journal = Journal.objects.first()
        consumer = ConsumerPubmed.objects.first()

        papers = Paper.objects.filter(journal=journal)
        self.assertTrue(papers.count() == 0)
        consumer.populate_journal(journal.id)
        papers = Paper.objects.filter(journal=journal)
        self.assertEqual(papers.count(), 10)


class ConsumerElsevierTest(TestCase):

    def setUp(self):
        journal, _ = Journal.objects.get_or_create(id_issn='1053-8119')
        consumer, _ = ConsumerElsevier.objects.get_or_create(name='pub cons #1')
        journal = Journal.objects.first()
        consumer.day0 = 10
        consumer.add_journal(journal)
        consumer.activate_journal(journal)
        consumer.save()

        with open('paperstream/consumers/tests/elsevier_sample.json') as file:
            self.entries = json.load(file)

    def test_add_entry(self):
        consumer = ConsumerElsevier.objects.first()
        journal = Journal.objects.first()
        item = consumer.parser.parse(self.entries[0])
        consumer.add_entry(item, journal)
        self.assertEqual(Paper.objects.all().count(), 1)

    @patch('paperstream.consumers.models.ConsumerElsevier.consume_journal')
    def test_consume_update_journal(self, mock_consume_journal):
        # setup mock object
        mock_consume_journal.return_value = self.entries, True

        journal = Journal.objects.first()
        consumer = ConsumerElsevier.objects.first()
        consumer.populate_journal(journal.pk)

        self.assertIsNotNone(consumer.consumerjournal_set.first().last_date_cons)
        self.assertTrue(
            consumer.consumerjournal_set.first().last_number_papers_recorded > 0)
        self.assertTrue(
            consumer.consumerjournal_set.first().last_number_papers_fetched > 0)

    @patch('paperstream.consumers.models.ConsumerElsevier.consume_journal')
    def test_populate_journal(self, mock_consume_journal):
        # setup mock object
        mock_consume_journal.return_value = self.entries, True

        journal = Journal.objects.first()
        consumer = ConsumerElsevier.objects.first()

        papers = Paper.objects.filter(journal=journal)
        self.assertTrue(papers.count() == 0)
        consumer.populate_journal(journal.pk)
        papers = Paper.objects.filter(journal=journal)
        self.assertEqual(papers.count(), 10)


class ConsumerArxivTest(TestCase):
# TODO: Implement
    pass