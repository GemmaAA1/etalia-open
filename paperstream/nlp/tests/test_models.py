import os
import tempfile
import shutil
import numpy as np

from django.core.exceptions import ValidationError
from django.conf import settings

from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES

from ..models import LSH, Model, TextField, JournalVectors, PaperVectors, \
    PaperNeighbors
from ..constants import FIELDS_FOR_MODEL
from ..exceptions import InvalidState

from .base import NLPTestCase, NLPDataTestCase


class ModelTest(NLPTestCase):

    def test_config_folders_have_been_created(self):
        self.assertTrue(os.path.isdir(settings.NLP_DOC2VEC_PATH))
        self.assertTrue(os.path.isdir(settings.NLP_DATA_PATH))
        self.assertTrue(os.path.isdir(settings.NLP_LSH_PATH))

    def test_model_can_be_instantiated(self):
        model = Model(name='test')
        model.full_clean()

    def test_model_print_name(self):
        model = Model(name='test')
        self.assertEqual(model.__str__(), 'test')

    def test_model_is_not_active_at_default(self):
        model = Model(name='test')
        self.assertFalse(model.is_active)

    def test_model_must_have_a_name(self):
        model = Model()
        with self.assertRaises(ValidationError):
            model.full_clean()

    def test_model_can_be_created(self):
        Model.objects.create(name='test')

    def test_models_cannot_have_same_name(self):
        Model.objects.create(name='test')
        with self.assertRaises(ValidationError):
            model = Model(name='test')
            model.full_clean()

    def test_doc2vec_model_can_be_saved(self):
        model = Model(name='test')
        model.save()
        path_to_be_save_to = os.path.join(
            settings.NLP_DOC2VEC_PATH,
            '{0}.mod'.format(model.name))
        self.assertTrue(os.path.isfile(path_to_be_save_to))

    def test_model_cannot_have_empty_text_field(self):
        with self.assertRaises(ValidationError):
            Model.objects.create(name='test', text_fields=[''])

    def test_model_has_default_text_fields(self):
        model = Model.objects.create(name='test')
        text_fields = [field for field in model.text_fields.all()]
        default = [field[0] for field in FIELDS_FOR_MODEL]
        self.assertTrue(set(default), set(text_fields))

    def test_model_can_have_different_fields(self):
        model = Model.objects.create(name='test', text_fields=['title'])
        model.full_clean()

    def test_model_can_have_different_size(self):
        model = Model(name='test', size=20)
        self.assertEqual(model.size, 20)

    def test_model_size_cannot_be_larger_than_NLP_MAX_VECTOR_SIZE(self):
        model = Model(name='test', size=settings.NLP_MAX_VECTOR_SIZE+1)
        with self.assertRaises(ValidationError):
            model.full_clean()


class ModelDumpTest(NLPDataTestCase):

    def test_model_can_dump_paper(self):
        self.model.dump(self.papers)

    def test_model_dumps_in_right_directory(self):
        self.model.dump(self.papers)
        self.assertTrue(os.path.isfile(
            os.path.join(settings.NLP_DATA_PATH,
                         '000000.txt')))

    def test_dump_file_start_with_paper_journal_pks(self):
        self.model.dump(self.papers)
        file = os.path.join(settings.NLP_DATA_PATH, '000000.txt')
        tag = '{pk}, j_{j_pk}:'.format(pk=self.papers[0].pk,
                                       j_pk=self.papers[0].journal.pk)
        with open(file, 'r') as f:
            line = f.readline()
            print(line[:len(tag)])
            print(tag)
            self.assertTrue(line[:len(tag)] == tag)

    def test_dump_file_splits_in_chunk(self):
        self.model.dump(self.papers)
        self.assertTrue(
            len([name for name in os.listdir(settings.NLP_DATA_PATH)]), 2)

    def test_dump_in_a_different_folder(self):
        data_path = tempfile.mkdtemp()
        self.model.dump(self.papers, data_path=data_path)
        self.assertTrue(os.path.isfile(os.path.join(data_path, '000000.txt')))
        # removing all file in there
        shutil.rmtree(data_path)

    def test_dump_file_is_trusted_false_discarded(self):
        self.model.dump([self.paper4])
        docs = self.model.load_documents()
        for count, _ in enumerate(docs):
            continue
        self.assertFalse('count' in locals())

    def test_dump_file_abstract_empty_discarded(self):
        self.model.dump([self.paper5])
        docs = self.model.load_documents()
        for count, _ in enumerate(docs):
            continue
        self.assertFalse('count' in locals())

    def test_dump_one_single_paper(self):
        self.model.dump([self.paper3])
        docs = self.model.load_documents()
        for count, _ in enumerate(docs):
            continue
        self.assertTrue('count' in locals())


class TestModelTrainingStack(NLPDataTestCase):

    def setUp(self):
        super(TestModelTrainingStack, self).setUp()
        self.model.dump(self.papers)

    def test_can_load_documents(self):
        self.model.load_documents()

    def test_can_build_phraser(self):
        docs = self.model.load_documents()
        self.model.build_phraser(docs, min_count=2)

    def test_can_build_vocab(self):
        docs = self.model.load_documents()
        self.model.build_vocab(docs)

    def test_can_train_save_and_be_set_active(self):
        docs = self.model.load_documents()
        self.model.build_vocab(docs)
        self.model.train(docs, passes=1, shuffle_=True)
        self.model.save()
        self.model.activate()
        self.assertTrue(self.model.is_active)

    def test_model_can_be_load_with_doc2vec(self):
        docs = self.model.load_documents()
        self.model.build_vocab_and_train()
        model2 = Model.objects.load(name='test')
        self.assertTrue(self.model == model2)

    def test_model_can_be_delete_with_doc2vec(self):
        docs = self.model.load_documents()
        self.model.build_vocab_and_train()
        self.model.save()
        self.model.delete()
        path_to_be_save_to = os.path.join(
            settings.NLP_DOC2VEC_PATH,
            '{0}.mod'.format(self.model.name))
        self.assertFalse(os.path.isfile(path_to_be_save_to))
        self.assertEqual(Model.objects.count(), 0)


class TestModelSavingStack(NLPDataTestCase):

    def setUp(self):
        super(TestModelSavingStack, self).setUp()
        self.model.dump(self.papers)
        self.model.build_vocab_and_train()

    def test_can_save_journal_vector(self):
        self.model.save_journal_vec_from_bulk()
        jv = JournalVectors.objects.get(journal=self.journal, model=self.model)
        self.assertIsNotNone(jv.vector)

    def test_can_save_paper_vector(self):
        self.model.save_paper_vec_from_bulk()
        pvs = PaperVectors.objects.all()
        for pv in pvs:
            self.assertIsNotNone(pv.vector)

    def test_can_infer_new_paper(self):
        self.model.infer_paper(paper_pk=self.paper.pk)
        pv = PaperVectors.objects.get(paper=self.paper, model=self.model)
        self.assertIsNotNone(pv.vector)

    # def test_can_infer_new_paper_with_seed(self):
    #     self.model.infer_paper(paper_pk=self.paper.pk, seed=True)
    #     pv = PaperVectors.objects.get(paper=self.paper, model=self.model)
    #     self.assertIsNotNone(pv.vector)

    def test_can_infer_multiple_new_paper(self):
        pks = [self.paper.pk, self.paper2.pk]
        self.model.infer_papers(pks)
        pv = PaperVectors.objects.get(paper=self.paper, model=self.model)
        self.assertIsNotNone(pv.vector)
        pv = PaperVectors.objects.get(paper=self.paper2, model=self.model)
        self.assertIsNotNone(pv.vector)

    def test_papervector_get_vector_length_match_model_size(self):
        self.model.infer_paper(paper_pk=self.paper.pk)
        pv = PaperVectors.objects.get(paper=self.paper, model=self.model)
        self.assertEqual(len(pv.get_vector()), self.model.size)


class TextFieldTest(NLPTestCase):

    def test_textfield_can_be_created(self):
        textfield = TextField(text_field=FIELDS_FOR_MODEL[0][0])
        textfield.full_clean()

    def test_textfield_must_be_among_choices(self):
        textfield = TextField(text_field='rqwerewqpmdspo')
        with self.assertRaises(ValidationError):
            textfield.full_clean()

    def test_textfield_cannot_be_empty(self):
        textfield = TextField(text_field='')
        with self.assertRaises(ValidationError):
            textfield.full_clean()

    def test_text_field_print_name(self):
        textfield = TextField(text_field=FIELDS_FOR_MODEL[0][0])
        self.assertEqual(textfield.__str__(), FIELDS_FOR_MODEL[0][0])


class PaperVectorTest(NLPDataTestCase):

    def test_papervector_cannot_be_created_if_model_not_active(self):
        self.model.is_active = False
        self.model.save_db_only()
        with self.assertRaises(ValidationError):
            PaperVectors.objects.create(paper=self.paper, model=self.model)

    def test_papervector_can_be_created_if_model_active(self):
        PaperVectors.objects.create(paper=self.paper, model=self.model)

    def test_papervector_must_be_unique_together_model_paper_pk(self):
        PaperVectors.objects.create(paper=self.paper, model=self.model)
        pv = PaperVectors(paper=self.paper, model=self.model)
        with self.assertRaises(ValidationError):
            pv.full_clean()

    def test_papervector_can_get_and_save_vector(self):
        pv = PaperVectors.objects.create(paper=self.paper, model=self.model)

        vector = np.random.randn(self.model.size)
        pv.set_vector(vector)

        pv2 = PaperVectors.objects.get(pk=pv.pk)

        vec_diff = np.sum(np.array(pv.get_vector()) -
                          np.array(pv2.get_vector()))

        self.assertAlmostEqual(vec_diff, 0)

class JournalVectorTest(NLPDataTestCase):

    def test_journalvector_cannot_be_created_if_model_NOT_active(self):
        self.model.deactivate()
        with self.assertRaises(ValidationError):
            JournalVectors.objects.create(journal=self.journal, model=self.model)

    def test_journalvector_can_be_created_if_model_active(self):
        self.model.activate()
        JournalVectors.objects.create(journal=self.journal, model=self.model)

    def test_journalvector_can_get_and_save_vector(self):
        jv = JournalVectors.objects.create(journal=self.journal,
                                           model=self.model)
        vector = np.random.randn(self.model.size)
        jv.set_vector(vector)

        jv2 = JournalVectors.objects.get(pk=jv.pk)

        vec_diff = np.sum(np.array(jv.get_vector()) -
                          np.array(jv2.get_vector()))

        self.assertAlmostEqual(vec_diff, 0)

    def test_journalvector_must_be_unique_together_model_paper_pk(self):
        JournalVectors.objects.create(journal=self.journal, model=self.model)
        jv = JournalVectors(journal=self.journal, model=self.model)
        with self.assertRaises(ValidationError):
            jv.full_clean()


class LSHModelTest(NLPDataTestCase):

    def setUp(self):
        super(LSHModelTest, self).setUp()
        self.pv = PaperVectors.objects.create(model=self.model,
                                              paper=self.paper)
        self.pv2 = PaperVectors.objects.create(model=self.model,
                                               paper=self.paper2)
        self.pv.set_vector(np.zeros(self.model.size))
        self.pv2.set_vector(np.zeros(self.model.size))

    def test_lsh_can_be_instantiated(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH.objects.create(model=self.model, time_lapse=time_lapse)

    def test_different_lsh_can_link_to_same_model(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH.objects.create(model=self.model, time_lapse=time_lapse)
        lsh = LSH.objects.create(model=self.model, time_lapse=-1)

    def test_lsh_default_state_is_none(self):
        model = Model.objects.first()
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=model, time_lapse=time_lapse)
        self.assertEqual(lsh.state, 'NON')

    def test_lsh_cannot_have_same_model_and_time_lapse(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        LSH.objects.create(model=self.model, time_lapse=time_lapse)
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        with self.assertRaises(ValidationError):
            lsh.full_clean()

    def test_lsh_time_lapse_raise_when_time_lapse_NOT_in_choices(self):
        lsh = LSH(model=self.model, time_lapse=1341241234312)
        with self.assertRaises(ValidationError):
            lsh.full_clean()

    def test_lsh_can_be_saved(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.save()
        save_path = os.path.join(settings.NLP_LSH_PATH,
                                 '{model_name}_{time_lapse}.lsh'
                                 .format(model_name=lsh.model.name,
                                         time_lapse=lsh.time_lapse))
        self.assertTrue(os.path.isfile(save_path))

    def test_lsh_can_be_loaded(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.save()
        lsh2 = LSH.objects.load(model=self.model, time_lapse=time_lapse)
        self.assertTrue(lsh, lsh2)

    def test_lsh_can_be_deleted(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.save()
        save_path = os.path.join(settings.NLP_LSH_PATH,
                                 '{model_name}_{time_lapse}.lsh'
                                 .format(model_name=lsh.model.name,
                                         time_lapse=lsh.time_lapse))
        lsh.delete()
        self.assertFalse(os.path.isfile(save_path))

    def test_lsh_can_change_state(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.set_state('BUS')
        self.assertEqual(lsh.state, 'BUS')

    def test_lsh_cannot_save_if_status_is_busy(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.set_state('BUS')
        with self.assertRaises(InvalidState):
            lsh.save()

    def test_lsh_get_data(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        x_data, pv_pks, new_pks = lsh.get_data()
        self.assertEqual(set(pv_pks), set([self.pv.pk, self.pv2.pk]))
        self.assertEqual(set(new_pks), set([self.paper.pk, self.paper2.pk]))
        self.assertTrue((x_data[0, :] == 0).all())

    def test_lsh_update_is_in_full_lsh_flag_if_full(self):
        time_lapse = -1
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        x_data, pv_pks, new_pks = lsh.get_data()
        lsh.update_if_full_lsh_flag(pv_pks)
        pv = PaperVectors.objects.get(pk=pv_pks[0])
        self.assertTrue(pv.is_in_full_lsh)

    def test_lsh_DOESNOT_update_is_in_full_lsh_flag_if_lsh_NOT_full(self):
        time_lapse = 61
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        x_data, pv_pks, new_pks = lsh.get_data()
        lsh.update_if_full_lsh_flag(pv_pks)
        pv = PaperVectors.objects.get(pk=pv_pks[0])
        self.assertFalse(pv.is_in_full_lsh)

    def test_lsh_can_run_update(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.update()

    def test_lsh_can_run_partial_update(self):
        time_lapse = -1
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.update(partial=True)

    def test_lsh_canNOT_run_partial_update_on_lsh_with_time_lapse(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        with self.assertRaises(ValueError):
            lsh.update(partial=True)

    def test_lsh_can_run_full_update(self):
        time_lapse = NLP_TIME_LAPSE_CHOICES[0][0]
        lsh = LSH(model=self.model, time_lapse=time_lapse)
        lsh.full_update()

    def test_lsh_can_be_created(self):
        LSH.objects.create(model=self.model,
                           time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])

    def test_lsh_can_search_for_kneighbors(self):
        lsh = LSH(model=self.model, time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        lsh.set_state('IDL')
        vec = np.zeros(lsh.model.size)
        indices = lsh.k_neighbors(vec)
        self.assertTrue(len(indices) > 0)

    def test_lsh_kneighbors_raises_if_state_not_idle(self):
        lsh = LSH(model=self.model, time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        vec = np.zeros(lsh.model.size)
        with self.assertRaises(InvalidState):
            lsh.k_neighbors(vec)

    def test_lsh_kneighbors_raises_with_wrong_vec_size(self):
        lsh = LSH(model=self.model, time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        lsh.set_state('IDL')
        vec = np.zeros(lsh.model.size+1)
        with self.assertRaises(AssertionError):
            lsh.k_neighbors(vec)

    def test_lsh_kneighbors_can_take_2d_array_in_input(self):
        lsh = LSH(model=self.model, time_lapse=NLP_TIME_LAPSE_CHOICES[0][0])
        lsh.set_state('IDL')
        mat = np.random.randn(10, lsh.model.size)
        indices = lsh.k_neighbors(mat, n_neighbors=2)
        self.assertEqual(indices.shape, (10, 2))


class PaperNeighborsTest(NLPDataTestCase):

    def setUp(self):
        super(PaperNeighborsTest, self).setUp()
        self.lsh = LSH(model=self.model, time_lapse=-1)
        self.lsh.save_db_only()

    def test_paperneighbors_can_be_instantiated(self):
        pn = PaperNeighbors(lsh=self.lsh, paper=self.paper)

    def test_paperneighbors_can_NOT_have_same_lsh_and_paper(self):
        PaperNeighbors.objects.create(lsh=self.lsh, paper=self.paper)
        pn = PaperNeighbors(lsh=self.lsh, paper=self.paper)
        with self.assertRaises(ValidationError):
            pn.full_clean()

    def test_paperneighbors_can_set_neighbors(self):
        pn = PaperNeighbors.objects.create(lsh=self.lsh, paper=self.paper)
        neigh = [1, 2, 3, 4]
        pn.set_neighbors(neigh)
        pn2 = PaperNeighbors.objects.get(pk=pn.id)
        self.assertTrue(pn2.get_neighbors() == pn.get_neighbors())


class LSHModelStackTest(NLPDataTestCase):

    def setUp(self):
        super(LSHModelStackTest, self).setUp()

    def test_can_build_all_lshs(self):
        self.model.dump(self.papers.all())
        self.model.build_vocab_and_train()
        self.model.save_journal_vec_from_bulk()
        self.model.save_paper_vec_from_bulk()
        self.model.build_lshs()

    def test_can_propagate(self):
        self.model.dump(self.papers.all())
        self.model.build_vocab_and_train()
        self.model.propagate()

    def test_can_build_full_stack(self):
        model2 = Model.objects.create(name='test2', size=32)
        model2.dump(self.papers.all())
        model2.build_vocab_and_train()
        model2.propagate()
