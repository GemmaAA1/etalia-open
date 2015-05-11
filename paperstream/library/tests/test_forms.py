from django.test import TestCase
from ..forms import JournalForm, PublisherForm
from ..models import Journal, Publisher


class PublisherFormTest(TestCase):

    def test_form_validation_for_blank_name(self):
        form = PublisherForm(data={'name': ''})
        self.assertFalse(form.is_valid())


class JournalFormTest(TestCase):

    def test_form_validation_can_save(self):
        publisher = Publisher.objects.create(name='Publisher #1')
        form = JournalForm(data={'title': 'Journal #1',
                                 'publisher': publisher.name})
        self.assertTrue(form.is_valid())

    def test_form_validation_for_blank_title(self):
        publisher = Publisher.objects.create(name='Publisher #1')
        form = JournalForm(data={'title': '', 'publisher': publisher.name})
        self.assertFalse(form.is_valid())

    def test_form_save_handles_saving_to_a_list(self):
        list_ = Publisher.objects.create(name)
        form = ItemForm(data={'text': 'do me'})
        new_item = form.save(for_list=list_)
        self.assertEqual(new_item, Item.objects.first())
        self.assertEqual(new_item.text, 'do me')
        self.assertEqual(new_item.list, list_)