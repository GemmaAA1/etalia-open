import numpy as np
from django.utils import timezone
from config.celery import celery_app as app
from core.constants import NLP_TIME_LAPSE_CHOICES
from library.models import Paper

from ..tasks import EmbedPaperTask, LSHTask, embed_all_models_and_find_neighbors, \
    embed_all_models, register_all_models_and_lshs_tasks
from ..models import Model, LSH
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


class LSHTaskClassTest(NLPDataExtendedTestCase):

    def setUp(self):
        super(LSHTaskClassTest, self).setUp()
        self.lsh = LSH.objects.get(model=self.model,
                                   time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])

    def test_lsh_task_can_be_instantiated(self):
        lsht = LSHTask(model_name=self.model.name,
                       time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        self.assertIsNone(lsht._lsh)

    def test_lsh_task_load_lsh_when_first_access(self):
        lsht = LSHTask(model_name=self.model.name,
                       time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        _ = lsht.lsh.state
        self.assertIsInstance(lsht._lsh, LSH)

    def test_lsh_task_model_name_must_exist(self):
        with self.assertRaises(ValueError):
            LSHTask(model_name='truc', time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])

    def test_lsh_task_time_lapse_must_be_among_choices(self):
        with self.assertRaises(ValueError):
            LSHTask(model_name='test', time_lapse=12312143431)

    def test_lsh_task_reload_lsh_if_modified(self):
        lsht = LSHTask(model_name=self.model.name,
                       time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        # load model
        state0 = lsht.lsh.state
        # modified lsh
        self.lsh.set_state('BUS')
        # test reload works
        self.assertEqual(self.lsh.state, lsht.lsh.state)
        # Test Not reload if model modification not saved
        self.lsh.state = 'IDL'
        self.assertNotEqual(self.lsh.state, lsht.lsh.state)
        # restore lsh state to idle
        self.lsh.set_state('IDL')


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
        register_all_models_and_lshs_tasks()

    def test_can_run_embed_paper_task(self):
        embed_task = app.tasks['nlp.tasks.embed_paper_{model_name}'
            .format(model_name=self.model.name)]
        self.assertIsNone(self.new_paper.vectors.first())
        rst = embed_task.delay(self.new_paper.pk)
        self.assertTrue(rst.successful())
        self.assertEqual(self.new_paper.vectors.count(), 1)
        self.assertIsNotNone(self.new_paper.vectors.first().get_vector())

    def test_can_run_all_embeddings_paper(self):
        embed_all_models(self.new_paper.pk)
        print(Model.objects.all().values_list('name', flat=True))
        self.assertEqual(self.new_paper.vectors.count(), 2)

    def test_vectors_have_correct_size(self):
        embed_all_models(self.new_paper.pk)
        vector_len = set([len(vec.get_vector())
                          for vec in self.new_paper.vectors.all()])
        size_len = set([model.size for model in Model.objects.all()])
        self.assertEqual(vector_len, size_len)


class LSHTaskTest(NLPDataExtendedTestCase):

    def setUp(self):
        super(LSHTaskTest, self).setUp()
        self.new_paper = Paper.objects.create(
            title='A New paper',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            date_ep=timezone.now().date())
        # train new model
        self.model2 = Model.objects.create(name='test2', size=32)
        self.model2.dump(self.papers.all())
        self.model2.build_vocab_and_train()
        self.model2.propagate()
        register_all_models_and_lshs_tasks()

    def test_update_task(self):
        lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_{time_lapse}'.format(
            model_name=self.model.name, time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])]
        res = lsh_task.delay('update')
        self.assertTrue(res.successful())

    def test_update_partial_task(self):
        lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_-1'.format(
            model_name=self.model.name)]
        res = lsh_task.delay('update')
        self.assertTrue(res.successful())

    def test_full_update_task(self):
        lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_-1'.format(
            model_name=self.model.name)]
        res = lsh_task.delay('full_update')
        self.assertTrue(res.successful())

    def test_populate_neighbors_task(self):
        lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_-1'.format(
            model_name=self.model.name)]
        res = lsh_task.delay(self.paper.pk, 'populate_neighbors')
        self.assertTrue(res.successful())

    def test_k_neighbors_task(self):
        lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_-1'.format(
            model_name=self.model.name)]
        pks = self.papers.values_list('pk', flat=True)
        model_pk = self.model.pk
        res = lsh_task.delay('k_neighbors_pks', seed_pks=pks,
                             model_pk=model_pk, k=4)
        self.assertTrue(res.successful())

    def test_embed_all_models_and_find_neighbors(self):
        embed_all_models_and_find_neighbors(self.new_paper.pk)
        self.assertEqual(self.new_paper.vectors.count(), 2)
        self.assertEqual(self.new_paper.neighbors.count(),
                         2 * len(NLP_TIME_LAPSE_CHOICES))
