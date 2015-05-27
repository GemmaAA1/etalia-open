from social.backends.oauth import BaseOAuth2
from social.backends.mendeley import MendeleyMixin


class CustomMendeleyOAuth2(MendeleyMixin, BaseOAuth2):

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
            out['affiliation'] = {}
            first_aff = response.get('employment')[0]
            first_aff_details = first_aff.get('institution_details', '')
            out['affiliation']['department'] = first_aff.get('department', '')
            out['affiliation']['institution'] = first_aff.get('institution', '')
            out['affiliation']['city'] = first_aff_details.get('city', '')
            out['affiliation']['state'] = first_aff_details.get('state', '')
            out['affiliation']['country'] = first_aff_details.get('country', '')
        return out
