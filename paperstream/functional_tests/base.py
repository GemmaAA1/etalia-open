import names
import random
# import factory
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from stdnum import issn
from library.models import Paper, Journal, Author, AuthorPosition
from library.forms import JournalForm, PaperForm, AuthorForm


class FunctionalTest(StaticLiveServerTestCase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.against_staging = False
        cls.server_url = cls.live_server_url

    @classmethod
    def tearDownClass(cls):
        if not cls.against_staging:
            super().tearDownClass()

    def setUp(self):
        self.browser = webdriver.Firefox()
        self.browser.implicitly_wait(3)

    def tearDown(self):
        self.browser.quit()

    def create_pre_load_library_with_journals_and_papers(self, *args):
        # if self.against_staging:
        #     pre_load_library_on_server(self.server_host)
        # else:
        #     pre_load_library()
        self.pre_load_library(*args)

    def pre_load_library(self, nb_journals=10, nb_authors=50, nb_papers=30):


        # Create a bench of journals
        journal_errors = []
        for n in range(nb_journals):
            digit7 = random.randint(0, 1e7-1)
            form = JournalForm(
                {'title': 'Journal title #{0:02d}'.format(n+1),
                 'id_issn': '{0:07d}{1}'.format(
                    digit7,
                    issn.calc_check_digit('{0:07d}'.format(digit7)))})
            # Check to see if the row data is valid.
            if form.is_valid():
                form.save()
            else:
                journal_errors.append(form.errors)


        # Create a bench of authors
        authors = []
        while len(authors) < nb_authors:
            authors.append(Author.objects.create(
                first_name=names.get_first_name(),
                last_name=names.get_last_name())
            )

        # Create a bench of papers with random authors
        papers = []
        n = 0
        while n < nb_papers:
            n += 1
            # build paper
            paper = Paper.objects.create(
                title='Title for paper #{0:02d}'.format(n),
                id_doi='doi.{0:03d}'.format(n),
                journal=journals[random.randint(0, nb_journals - 1)])
            # add authors
            nb_authors_paper = random.randint(0, 5)
            for a in range(nb_authors_paper):
                AuthorPosition.objects.create(
                    paper=paper,
                    author=authors[random.randint(0, nb_authors - 1)],
                    position=a)



## Cant make factory_boy work with subfactory ref... :(

# class JournalFactory(factory.Factory):
#     class Meta:
#         model = Journal
#
#     title = factory.Sequence(lambda n: 'Journal title #{0:02d}'.format(n))
#     issn = factory.Sequence(lambda n: '{0:07d}{1}'.format(
#         n,
#         issn.calc_check_digit('{0:07d}'.format(n))))
#     ext_id = ''
#
#
# class AuthorFactory(factory.Factory):
#     class Meta:
#         model = Author
#
#     first_name = factory.Sequence(lambda n: '{0}'.format(names.get_first_name()))
#     last_name = factory.Sequence(lambda n: '{0}'.format(names.get_last_name()))
#
#
# class PaperFactory(factory.Factory):
#     class Meta:
#         model = Paper
#
#     title = factory.Sequence(lambda n: 'Title of paper #{0:02d}'.format(n))
#     journal = factory.SubFactory(JournalFactory)
#     doi = factory.Sequence(lambda n: 'doi.{0:03d}'.format(n))
#     authors = None
#
#
# class AuthorPositionFactory(factory.Factory):
#     class Meta:
#         model = AuthorPosition
#
#     paper = factory.SubFactory(PaperFactory)
#     author = factory.SubFactory(AuthorFactory)
#     position = 0