from social.backends.oauth import BaseOAuth2
from social.backends.mendeley import MendeleyMixin
from mendeley import Mendeley
from mendeley.auth import MendeleySession, \
    MendeleyAuthorizationCodeTokenRefresher
from mendeley.exception import MendeleyApiException

from .BaseMixin import BackendLibMixin
from .parsers import ParserMendeley

class CustomMendeleyOAuth2(MendeleyMixin, BackendLibMixin, BaseOAuth2):

    # library specific
    _type = 'MEN'
    parser = ParserMendeley()

    # base
    name = 'custom-mendeley-oauth2'
    REQUIRES_EMAIL_VALIDATION = True
    AUTHORIZATION_URL = 'https://api-oauth2.mendeley.com/oauth/authorize'
    ACCESS_TOKEN_URL = 'https://api-oauth2.mendeley.com/oauth/token'
    ACCESS_TOKEN_METHOD = 'POST'
    DEFAULT_SCOPE = ['all']
    REDIRECT_STATE = False
    EXTRA_DATA = MendeleyMixin.EXTRA_DATA + [
        ('refresh_token', 'refresh_token'),
        ('expires_in', 'expires_in'),
        ('expires_at', 'expires_at'),
        ('token_type', 'token_type'),
    ]



    def get_user_data(self, access_token, *args, **kwargs):
        """Loads user data from service"""
        return self.get_json(
            'https://api.mendeley.com/profiles/me/',
            headers={'Authorization': 'Bearer {0}'.format(access_token)}
        )

    def get_user_details(self, response):
        """Return user details from Mendeley account"""
        out = super(CustomMendeleyOAuth2, self).get_user_details(response)
        out['email'] = response.get('email', '')
        out['first_name'] = response.get('first_name', '')
        out['last_name'] = response.get('last_name', '')
        if response.get('employment', None):
            out['tmp_affiliation'] = {}
            first_aff = response.get('employment')[0]
            first_aff_details = first_aff.get('institution_details', '')
            out['tmp_affiliation']['department'] = first_aff.get('department', '')
            out['tmp_affiliation']['institution'] = first_aff.get('institution', '')
            out['tmp_affiliation']['city'] = first_aff_details.get('city', '')
            out['tmp_affiliation']['state'] = first_aff_details.get('state', '')
            out['tmp_affiliation']['country'] = first_aff_details.get('country', '')
        return out

    def get_session(self, social, user, *args, **kwargs):

        key, secret = self.get_key_and_secret()

        tokens = {'access_token': social.extra_data['access_token'],
                  'refresh_token': social.extra_data['refresh_token']}

        # start authorization flow
        mendeley = Mendeley(key, client_secret=secret,
                            redirect_uri=self.redirect_uri)

        auth = mendeley.start_authorization_code_flow()

        # start mendeley session
        mendeley_session = MendeleySession(
            mendeley, tokens, client=auth.client,
            refresher=MendeleyAuthorizationCodeTokenRefresher(auth))

        # try session
        try:
            name = mendeley_session.profiles.me.display_name
        except MendeleyApiException:
            # renew token
            new_token = self.refresh_token(social.extra_data['refresh_token'])
            social.extra_data['access_token'] = new_token['access_token']
            social.save()
            tokens['access_token'] = social.extra_data['access_token']
            # start mendeley session
            mendeley_session = MendeleySession(
                mendeley, tokens, client=auth.client,
                refresher=MendeleyAuthorizationCodeTokenRefresher(auth))

        return mendeley_session

    def update_lib(self, mendeley_session, user, *args, **kwargs):

        # Init
        new = True
        count = 0
        user.lib.set_lib_syncing()

        # retrieve list of documents per page
        page = mendeley_session.documents.list(
            page_size=self.CHUNK_SIZE,
            sort='created',
            order='desc',
            view='all')

        # loop through documents/page until end or find already registered doc
        while True:
            for item in page.items:
                entry = self.parser.parse(item)
                paper, journal = self.add_entry(entry)
                if paper:
                    new = self.associate_paper(paper, user, entry['user_info']) and new
                    if new:
                        count += 1
                    else:
                        break
                    if journal:
                        self.associate_journal(journal, user)
            if page.next_page:
                page = page.next_page
            else:
                break



        # update UserLib and Stats
        self.create_lib_stats(user, count)
        user.lib.set_lib_idle()
        return count





