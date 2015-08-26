import os
import shutil
from django.conf import settings
from django.test import TestCase
from django.utils import timezone

from paperstream.library.models import Journal, Paper

from ..models import Model
from ..tasks import register_all_models_and_lshs_tasks

class NLPTestCase(TestCase):

    def tearDown(self):
        # Removing test folders for NLP
        if os.path.isdir(settings.NLP_DOC2VEC_PATH):
            shutil.rmtree(settings.NLP_DOC2VEC_PATH)

        if os.path.isdir(settings.NLP_DATA_PATH):
            shutil.rmtree(settings.NLP_DATA_PATH)

        if os.path.isdir(settings.NLP_LSH_PATH):
            shutil.rmtree(settings.NLP_LSH_PATH)



class NLPDataTestCase(NLPTestCase):

    def setUp(self):
        self.journal = Journal.objects.create(title='Journal test')
        self.paper = Paper.objects.create(
            title='Bla bla.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=True)
        self.paper2 = Paper.objects.create(
            title='Blo blo.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=True)
        self.paper3 = Paper.objects.create(
            title='Bli bli.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=True)
        self.paper4 = Paper.objects.create(
            title='Bli bli.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=False)
        self.paper5 = Paper.objects.create(
            title='Bli bli.',
            abstract='',
            journal=self.journal,
            date_ep=timezone.now().date(),
            is_trusted=True)

        self.model = Model.objects.create(name='test',
                                          size=128)
        self.model.activate()
        self.model.save_db_only()
        self.papers = Paper.objects.all()


class NLPDataExtendedTestCase(NLPDataTestCase):

    def setUp(self):
        super(NLPDataExtendedTestCase, self).setUp()
        self.model.dump(self.papers.all())
        self.model.build_vocab_and_train()
        self.model.propagate()
        register_all_models_and_lshs_tasks()
