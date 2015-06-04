from .base import FunctionalTest


class LoginTest(FunctionalTest):

    def test_login_with_mendeley(self):

        # X lands on paperstream home page and notice mendeley login
        self.browser.get(self.server_url)
        self.browser.find_element_by_id('mendeley_login').click()

        # A Mendeley login box appears
        self.switch_to_new_window('Mendeley')

        # X login with it email and password

        # The mendeley window closes

        # X is redirected to a pre-filed form that need to be completed

        # X completes the form and sign-in

        # X is redirected to page that ask him to check his email to activate
        # his account

        # In the meantime the importation of his library has started

        # X clicks on the activated link in his email and is redirected to
        # a page telling him that is now activated

        # ...