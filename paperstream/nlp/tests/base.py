import os
import shutil
from django.conf import settings
from django.test import TestCase
from django.utils import timezone

from library.models import Journal, Paper
from ..models import LSH, Model, PaperVectors

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
            date_ep=timezone.now().date())
        self.paper2 = Paper.objects.create(
            title='Blo blo.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date())
        self.paper3 = Paper.objects.create(
            title='Bli bli.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date())
        self.model = Model.objects.create(name='test',
                                          size=128)
        self.model.activate()
        self.model.save_db_only()
        self.papers = Paper.objects.all()
