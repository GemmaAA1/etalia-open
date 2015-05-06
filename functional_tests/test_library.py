from .base import FunctionalTest
from library.models import Journal, Paper

class LibraryTest(FunctionalTest):

    def test_browse_library_anonymous(self):

        # Paper stream is already loaded with some papers and journals
        self.create_pre_load_library_with_journals_and_papers()

        # Anonymous X wants to check what the library has in stock.
        # He goes to the root library URL
        self.browser.get(self.server_url+'/library/')
        library_landing_url = self.browser.current_url

        # Th library landing page of the library shows some statistics
        # journals, papers, creators

        # Then he clicks on journal
        self.browser.find_element_by_link_text('Journal').click()

        # X is redirected to library/journals/
        self.assertEqual(self.browser.current_url,
                         self.server_url+'/library/journals/')

        # Where a list of journals in db is displayed. First one is titled
        # "Journal title #01"

        # X clicks on the first Journal
        # self.browser.find_element_by_link_text('Journal title #01').click()

        # X is redirected to /library/journals/<pk>/papers/
        # pk = Journal.objects.get(title='Journal title #01').id
        # self.assertEqual(self.browser.current_url,
        #                  '{0}/library/journals/{1}/papers'.format(
        #                      self.server_url, pk))

        # The page displays a list of paper for this journal
        # p = Paper.objects.filter(journal__pk=pk).first()


        # X clicks on the first paper of the list

        # X is redirected to the /library/paper/<pk> where attributes of the
        # paper are displayed
