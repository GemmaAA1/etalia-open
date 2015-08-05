from django.test import TestCase

from core.constants import TIME_LAPSE_CHOICES

from ..tasks import EmbedPaperTask, LSHTask, all_embeddings_and_neighbors, \
    all_embedings_paper
from ..models import Model, LSH
from .base import NLPDataTestCase, NLPDataExtendedTestCase

class EmbedPaperTaskTest(TestCase):

    def setUp(self):
        self.model = Model.objects.create(name='test')

    def test_embed_paper_task_can_be_instantiated(self):
        ept = EmbedPaperTask(model_name=self.model.name, bind=True)
        self.assertIsNone(ept._model)

    def test_embed_paper_task_can_load_model_when_first_access(self):
        ept = EmbedPaperTask(model_name=self.model.name, bind=True)
        _ = ept.model.name
        self.assertIsInstance(ept._model, Model)

    def test_embed_paper_task_must_have_a_model_name(self):
        with self.assertRaises(TypeError):
            EmbedPaperTask()

    def test_embed_paper_task_model_name_must_be_defined(self):
        with self.assertRaises(ValueError):
            EmbedPaperTask(model_name='truc')

    def test_embed_paper_task_reload_model_if_modified(self):
        ept = EmbedPaperTask(model_name=self.model.name, bind=True)
        # load model
        size0 = ept.model.size
        # modified model
        self.model.size = int(self.model.size / 2)
        self.model.save()
        # test reload works
        self.assertEqual(self.model.size, ept.model.size)
        # Test Not reload if model modification not saved
        self.model.size = int(self.model.size * 2)
        self.assertNotEqual(self.model.size, ept.model.size)


class LSHTaskTest(NLPDataExtendedTestCase):

    def setUp(self):
        super(LSHTaskTest, self).setUp()
        self.lsh = LSH.objects.create(model=self.model,
                                      time_lapse=TIME_LAPSE_CHOICES[0][0])

    def test_lsh_task_can_be_instantiated(self):
        lsht = LSHTask(model_name=self.model.name,
                       time_lapse=TIME_LAPSE_CHOICES[0][0])
        self.assertIsNone(lsht._lsh)

    def test_lsh_task_load_lsh_when_first_access(self):
        lsht = LSHTask(model_name=self.model.name,
                       time_lapse=TIME_LAPSE_CHOICES[0][0])
        _ = lsht.lsh.state
        self.assertIsInstance(lsht._lsh, LSH)

    def test_lsh_task_must_have_model_name_and_time_lapse(self):
        with self.assertRaises(TypeError):
            LSHTask()
        with self.assertRaises(TypeError):
            LSHTask(model_name=self.model.name)
        with self.assertRaises(TypeError):
            LSHTask(time_lapse=TIME_LAPSE_CHOICES[0][0])

    def test_lsh_task_model_name_must_exist(self):
        with self.assertRaises(ValueError):
            LSHTask(model_name='truc', time_lapse=TIME_LAPSE_CHOICES[0][0])

    def test_lsh_task_time_lapse_must_be_among_choices(self):
        with self.assertRaises(ValueError):
            LSHTask(model_name='test', time_lapse=12312143431)

