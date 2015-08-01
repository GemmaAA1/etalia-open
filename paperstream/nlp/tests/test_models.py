from django.test import TestCase
from django.core.exceptions import ValidationError
from ..models import LSH, Model
from library.models import Paper


class ModelModelTest(TestCase):

    def setUp(self):
        pass

    def model_can_be_create(self):
        Model.objects.create(name='test')

    def model_is_not_active_at_default(self):
        model = Model(name='test')
        self.assertFalse(model.is_active)

    def model_must_have_a_name(self):
        with self.assertRaises(ValidationError):
            model = Model()
            model.full_clean()

    def models_cannot_have_same_name(self):
        Model.objects.create(name='test')
        with self.assertRaises(ValidationError):
            model = Model(name='test')
            model.full_clean()



class lSHModelTest(TestCase):

    def setUp(self):
        Model.objects.create(name='test')
        Paper.objects.create(title='On the road again')

    def test_lsh_can_be_create(self):
        pass