import names
import random
# import factory
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from stdnum import issn
from library.models import Paper, Journal, Author, AuthorPaper


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
        journals = []
        for n in range(nb_journals):
            digit7 = random.randint(0, 1e7-1)
            journals.append(Journal.objects.create(
                title='Journal title #{0:02d}'.format(n+1),
                id_issn='{0:07d}{1}'.format(digit7, issn.calc_check_digit(
                    '{0:07d}'.format(digit7))))
            )

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
                AuthorPaper.objects.create(
                    paper=paper,
                    author=authors[random.randint(0, nb_authors - 1)],
                    position=a)