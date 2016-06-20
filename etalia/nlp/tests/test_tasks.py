# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.utils import timezone

from config.celery import celery_app as app
from etalia.core.constants import NLP_TIME_LAPSE_CHOICES
from etalia.library.models import Paper
from library.tasks import embed_paper

from ..tasks import EmbedPaperTask, MostSimilarTask, register_nlp_tasks

from ..models import Model, PaperEngine
from .base import NLPDataExtendedTestCase, NLPDataTestCase


class EmbedPaperTaskClassTest(NLPDataTestCase):

    def setUp(self):
        super(EmbedPaperTaskClassTest, self).setUp()
        self.model.dump(self.papers.all())
        self.model.build_vocab_and_train()
        self.model.propagate()

    def test_embed_paper_task_can_be_instantiated(self):
        ept = EmbedPaperTask(model_name=self.model.name, bind=True)
        self.assertIsNone(ept._model)

    def test_embed_paper_task_can_load_model_when_first_access(self):
        ept = EmbedPaperTask(model_name=self.model.name, bind=True)
        _ = ept.model.name
        self.assertIsInstance(ept._model, Model)

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


class EmbedTaskTest(NLPDataExtendedTestCase):

    def setUp(self):
        super(EmbedTaskTest, self).setUp()
        self.new_paper = Paper.objects.create(
            title='A New paper',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            date_ep=timezone.now().date())
        # train new model
        self.model2 = Model.objects.create(name='test2', size=32)
        self.model2.dump(self.papers)
        self.model2.build_vocab_and_train()
        register_nlp_tasks()

    def test_can_run_all_embeddings_paper(self):
        embed_paper(self.new_paper.pk)
        print(Model.objects.all().values_list('name', flat=True))
        self.assertEqual(self.new_paper.vectors.count(), 2)

    def test_vectors_have_correct_size(self):
        embed_paper(self.new_paper.pk)
        vector_len = set([len(vec.get_vector())
                          for vec in self.new_paper.vectors.all()])
        size_len = set([model.size for model in Model.objects.all()])
        self.assertEqual(vector_len, size_len)
