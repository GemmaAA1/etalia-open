import json

from allauth.socialaccount.providers.oauth.client import OAuth
from allauth.socialaccount.providers.oauth.views import (OAuthAdapter,
                                                         OAuthLoginView,
                                                         OAuthCallbackView)
from .provider import ZoteroProvider



class ZoteroAPI(OAuth):

    def get_user_info(self):
        data = self._get_at_from_session()
        return data

class ZoteroOAuthAdapter(OAuthAdapter):
    provider_id = ZoteroProvider.id
    request_token_url = 'https://www.zotero.org/oauth/request'
    access_token_url = 'https://www.zotero.org/oauth/access'
    authorize_url = 'https://www.zotero.org/oauth/authorize'

    def complete_login(self, request, app, token):
        client = ZoteroAPI(request, app.client_id, app.secret,
                           self.request_token_url)
        extra_data = client.get_user_info()

        return self.get_provider().sociallogin_from_response(request,
                                                             extra_data)

oauth_login = OAuthLoginView.adapter_view(ZoteroOAuthAdapter)
oauth_callback = OAuthCallbackView.adapter_view(ZoteroOAuthAdapter)
