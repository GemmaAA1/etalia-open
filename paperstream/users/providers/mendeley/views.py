import requests
from allauth.socialaccount.providers.oauth2.views import (OAuth2Adapter,
                                                          OAuth2LoginView,
                                                          OAuth2CallbackView)
from .provider import MendeleyProvider


class MendeleyOAuth2Adapter(OAuth2Adapter):
    provider_id = MendeleyProvider.id
    access_token_url = 'https://api.mendeley.com/oauth/token'
    authorize_url = 'https://api.mendeley.com/oauth/authorize'
    profile_url = 'https://api.mendeley.com/profiles/me'

    def complete_login(self, request, app, token, **kwargs):
        resp = requests.get(self.profile_url,
                            params={'access_token': token.token})
        extra_data = resp.json()

        return self.get_provider().sociallogin_from_response(request,
                                                             extra_data)


oauth2_login = OAuth2LoginView.adapter_view(MendeleyOAuth2Adapter)
oauth2_callback = OAuth2CallbackView.adapter_view(MendeleyOAuth2Adapter)

