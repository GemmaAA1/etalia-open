from .base import FunctionalTest

MENDELEY_EMAIL = 'trucfortest@gmail.com'
MENDELEY_PASSWORD = 'qwerty123'

class AuthenticateTest(FunctionalTest):

    def test_signup_with_mendeley(self):

        # X lands on paperstream home page and notice mendeley login
        self.browser.get(self.server_url)
        self.browser.find_element_by_id('mendeley_login').click()

        # He is redirected to the mendeley auth page
        self.switch_to_new_window('Mendeley')
        # He enters its login/password
        self.browser.find_element_by_id('username').send_keys()
        self.browser.find_element_by_tag_name('button').click()


        # X is redirected to sign up profile page with prepopulated field

        # He presses continue and an email is send out

        # Its library starts to be uploaded in the backgroud

        # he clicks on the link

        # He is asked to fill up his affiliation (pre-populated)

        # He hits continue and is redirected to landing authenticated guys.

        # X is redirected to a pre-filed form that need to be completed

        # X completes the form and sign-in

        # X is redirected to page that ask him to check his email to activate
        # his account

        # In the meantime the importation of his library has started

        # X clicks on the activated link in his email and is redirected to
        # a page telling him that is now activated

        # ...


    def test_signup_with_zotero(self):

        # X lands on paperstream home page and notice mendeley login
        self.browser.get(self.server_url)
        self.browser.find_element_by_id('zotero_login').click()

        # He is redirected to the mendeley auth page
        self.switch_to_new_window('Mendeley')
        # He enters its login/password
        self.browser.find_element_by_id('authentication_email'
                                        ).send_keys(TEST_EMAIL)
        self.browser.find_element_by_tag_name('button').click()
