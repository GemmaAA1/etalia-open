from .base import FunctionalTest
from unittest import skip

MENDELEY_EMAIL = 'trucfortest@gmail.com'
MENDELEY_FIRST_NAME = 'Toto'
MENDELEY_LAST_NAME = 'Popo'
MENDELEY_PASSWORD = 'qwerty123'
ZOTERO_EMAIL = 'trucfortest@gmail.com'
ZOTERO_PASSWORD = 'qwerty123'

@skip
class AuthenticateWithZoteroTest(FunctionalTest):

    # def test_signup_with_zotero(self):
    #
    #     # X lands on paperstream home page and notice mendeley login
    #     self.browser.get(self.server_url)
    #     self.browser.find_element_by_id('zotero-link').click()
    #
    #     # X is redirected to the zotero login page
    #     # X enters its login/password and click log in
    #     self.browser.find_element_by_id('username').send_keys(ZOTERO_EMAIL)
    #     self.browser.find_element_by_id('password').send_keys(ZOTERO_PASSWORD)
    #     self.browser.find_element_by_id('login').click()
    #
    #     # X is redirected to private key page and X accepts
    #     self.browser.find_element_by_id('accept').click()
    #
    #     # X is redirected to the basic info page
    #     self.browser.find_element_by_id('email').send_keys(ZOTERO_EMAIL)
    #     self.browser.find_element_by_id('first_name').send_keys('Zotero')
    #     self.browser.find_element_by_id('last_name').send_keys('Zotero')
    #     self.browser.find_element_by_class_name('btn btn-primary btn-lg').click()




class AuthenticationWithMendeleyTest(FunctionalTest):

    def test_signup_with_mendeley(self):

        # X lands on paperstream home page and notice mendeley login
        self.browser.get(self.server_url)
        self.browser.find_element_by_id('mendeley-link').click()

        # It is redirected to the mendeley auth page
        self.switch_to_new_window('Mendeley')
        # It enters its login/password
        self.browser.find_element_by_id('username').send_keys(MENDELEY_EMAIL)
        self.browser.find_element_by_id('password').send_keys(MENDELEY_PASSWORD)
        self.browser.find_element_by_class_name('btn btn-primary').click()

        # It is redirected to user/info
        self.wait_for_element_with_id('id_first_name')
        current_url = self.browser.current_url()

        # check callback
        self.assertIn('user/info', current_url)
        # check mendeley is returning email / first_name / last_name
        first_name = self.browser.find_element_by_name('first_name').value
        last_name = self.browser.find_element_by_name('last_name').value
        email = self.browser.find_element_by_name('email').value
        self.assertEqual(first_name, MENDELEY_FIRST_NAME)
        self.assertEqual(last_name, MENDELEY_LAST_NAME)
        self.assertEqual(email, MENDELEY_EMAIL)
        # click continue
        self.browser.find_element_by_xpath(
            '//*[@id="modalPrimaryInfo"]/div/div/div[2]/form/div[4]/button')\
            .click()

        #

#
#         X is redirected to sign up profile page with prepopulated field
#
#         He presses continue and an email is send out
#
#         Its library starts to be uploaded in the backgroud
#
#         he clicks on the link
#
#         He is asked to fill up his affiliation (pre-populated)
#
#         He hits continue and is redirected to landing authenticated guys.
#
#         X is redirected to a pre-filed form that need to be completed
#
#         X completes the form and sign-in
#
#         X is redirected to page that ask him to check his email to activate
#         his account
#
#         In the meantime the importation of his library has started
#
#         X clicks on the activated link in his email and is redirected to
#         a page telling him that is now activated

        # ...