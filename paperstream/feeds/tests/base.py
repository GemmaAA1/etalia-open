from django.test import TestCase
from django.utils import timezone
import numpy as np

from library.models import Journal, Paper
from nlp.models import LSH, Model, PaperVectors, JournalVectors
from users.models import User, UserLib, UserLibPaper, UserSettings
from nlp.tasks import register_all_models_and_lshs_tasks

class UserFeedTestCase(TestCase):

    def tearDown(self):
        pass

    def setUp(self):
        self.journal = Journal.objects.create(title='Journal test')
        self.journal1 = Journal.objects.create(title='Journal test 2')
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
            journal=self.journal1,
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
        self.paper6 = Paper.objects.create(
            title='Bli bli.',
            abstract='',
            journal=self.journal1,
            date_ep=timezone.now().date(),
            is_trusted=True)

        self.model = Model.objects.create(name='test',
                                          size=128)
        self.model.is_active = True
        self.model.save_db_only()
        self.user = User(email='test@paperstream.com')
        self.user.save()
        UserSettings.objects.create(user=self.user)
        self.papers = Paper.objects.all()
        self.journals = Journal.objects.all()
        ul = UserLib(user=self.user)
        ul.save()
        # Papers 4, 5, 6 go in user lib
        UserLibPaper.objects.create(
            userlib=ul,
            paper=self.paper6,
            date_created=(timezone.now() - timezone.timedelta(days=300)).date())
        UserLibPaper.objects.create(
            userlib=ul,
            paper=self.paper5,
            date_created=(timezone.now() - timezone.timedelta(days=100)).date())
        UserLibPaper.objects.create(
            userlib=ul,
            paper=self.paper4,
            date_created=(timezone.now() - timezone.timedelta(days=10)).date())

        # add random vector to papers
        for paper in self.papers:
            pv = PaperVectors.objects.create(paper=paper, model=self.model)
            vec = np.random.randn(self.model.size)
            pv.set_vector(vec)

        # add random vector to journals
        for journal in self.journals:
            jv = JournalVectors.objects.create(journal=journal, model=self.model)
            vec = np.random.randn(self.model.size)
            jv.set_vector(vec)

        # register celery task for testing
        self.model.build_lshs()
        register_all_models_and_lshs_tasks()