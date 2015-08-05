from django.test import TestCase

from ..tasks import EmbedPaperTask, LSHTask, all_embeddings_and_neighbors, \
    all_embedings_paper

from ..models import Model

class EmbedTest(TestCase):

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
