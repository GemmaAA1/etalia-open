import time
from unittest import skip
from .base import FunctionalTest
from django.conf import settings
from library.models import Journal, Paper

@skip
class LibraryTest(FunctionalTest):

    def test_browse_library_anonymous(self):

        # Paper stream is already loaded with some papers and journals
        self.create_pre_load_library_with_journals_and_papers(
            settings.ITEMS_PER_PAGE * 2,
            50,
            settings.ITEMS_PER_PAGE * 10)

        # Anonymous X wants to check what the library has in stock.
        # He goes to the root library URL
        self.browser.get(self.server_url+'/library/')
        library_landing_url = self.browser.current_url

        # Th library landing page of the library shows some statistics
        # journals, papers, creators

        # Then he clicks on journal
        self.browser.find_element_by_link_text('Journals').click()

        # X is redirected to library/journals/
        self.assertEqual(self.browser.current_url,
                         self.server_url+'/library/journals/')

        # X checks next page
        self.browser.find_element_by_link_text('Next').click()

        # X comes back to previous page
        self.browser.find_element_by_link_text('Previous').click()

        # X clicks on the first Journal
        self.browser.find_element_by_id('journal1').click()

        # X is redirected to /library/journal/<pk>/
        pk = Journal.objects.get(title='Journal title #01').id
        self.assertEqual(self.browser.current_url,
                         '{0}/library/journal/{1}/'.format(
                             self.server_url, pk))

        # The page displays a list of paper for this journal
        p = Paper.objects.filter(journal__pk=pk).first()

        # X clicks on the first paper of the list
        self.browser.find_element_by_id('paper1').click()

        # X is redirected to the /library/paper/<pk> where attributes of the
        # paper are displayed
        pk = Paper.objects.get(title='Journal title #01').id
        self.assertEqual(self.browser.current_url,
                         '{0}/library/paper/{1}/'.format(
                             self.server_url, pk))
