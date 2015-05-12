from django.test import TestCase
from ..forms import JournalForm, PublisherForm
from ..models import Journal, Publisher


class PublisherFormTest(TestCase):

    def test_form_validation(self):
        form = PublisherForm(data={'name': 'Publisher Zorro'})
        self.assertTrue(form.is_valid())

    def test_form_validation_blank_name_invalid(self):
        form = PublisherForm(data={'name': ''})
        self.assertFalse(form.is_valid())


class JournalFormTest(TestCase):

    def test_form_validation(self):
        form = JournalForm(data={'title': ''})
        self.assertFalse(form.is_valid())

    def test_form_clean_title(self):
        form = JournalForm(data={'title': 'JOURNAL CAPITALS'})
        self.assertTrue(form.is_valid())
        self.assertEqual('Journal Capitals', form.cleaned_data['title'])

    def test_form_journal_has_publisher_and_save(self):
        publisher = Publisher.objects.create(name='Publisher #1')
        form = JournalForm(data={'title': 'Journal #1',
                                 'publisher': publisher.name})
        self.assertTrue(form.is_valid())
        form.save()
        publisher_saved = Journal.objects.first().publisher
        self.assertEqual(publisher_saved, publisher)

    def test_form_validation_blank_title_is_invalid(self):
        publisher = Publisher.objects.create(name='Publisher #1')
        form = JournalForm(data={'title': '', 'publisher': publisher.name})
        self.assertFalse(form.is_valid())
