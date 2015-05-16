from django.test import TestCase
from ..models import ConsumerPubmed, ConsumerJournal
from library.models import Journal
from django.core.exceptions import ValidationError


class ConsumerTest(TestCase):

    def setUp(self):
        Journal.objects.create(title='Journal of the consumers',
                               id_issn='0000-0019')

    @staticmethod
    def add_journal():
        consumer = ConsumerPubmed.objects.create()
        journal = Journal.objects.first()
        consumer.add_journal(journal)
        cj = consumer.add_journal(journal)
        return consumer, journal, cj

    def test_consumer_have_unique_name(self):
        ConsumerPubmed.objects.create(name='consumer #1')
        with self.assertRaises(ValidationError):
            ConsumerPubmed(name='consumer #1').full_clean()

    def test_can_add_and_retrieve_journal(self):
        consumer = ConsumerPubmed.objects.create()
        journal = Journal.objects.first()
        consumer.add_journal(journal)
        cj = ConsumerJournal.objects.all()[0]
        self.assertEqual(journal, cj.journal)
        self.assertEqual(consumer.pk, cj.consumer.pk)
        self.assertEqual(cj.status, 'inactive')

    def test_can_activate_journal(self):
        consumer, journal, cj = self.add_journal()
        consumer.activate_journal(journal)
        cj = consumer.consumerjournal_set.get(journal=journal)
        self.assertEqual(cj.status, 'idle')

    def test_can_deactivate_journal(self):
        consumer, journal, cj = self.add_journal()
        consumer.activate_journal(journal)
        consumer.deactivate_journal(journal)
        self.assertEqual(cj.status, 'inactive')

    def test_cannot_activate_if_error_status(self):
        consumer, journal, cj = self.add_journal()
        cj.status = 'error'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.activate_journal(journal)

    def test_cannot_activate_if_doing_status(self):
        consumer, journal, cj = self.add_journal()
        cj.status = 'consuming'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.activate_journal(journal)
        cj.status = 'in_queue'
        cj.save()
        with self.assertRaises(ValueError):
            consumer.activate_journal(journal)

    def test_cannot_deactivate_if_doing_stuff(self):
        consumer, journal, cj = self.add_journal()
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
        journal = Journal.objects.first()
        consumer.add_journal(journal)

    def test_consume_journal_update_consumerjournal(self):
        journal = Journal.objects.first()
        consumer = ConsumerPubmed.objects.first()
        consumer.consume(journal)
        self.assertIsNotNone(consumer.consumerjournal_set.first().last_date_cons)
        self.assertTrue(consumer.consumerjournal_set.first().last_number_papers > 0)

    def test_consume_journal_update_consumerjournalstat(self):
        journal = Journal.objects.first()
        consumer = ConsumerPubmed.objects.first()
        consumer.consume(journal)
        consumer.consume(journal)
        cj = consumer.consumerjournal_set.first()
        cjs = cj.stats_set.all()
        self.assertEqual(cjs.count(), 2)
