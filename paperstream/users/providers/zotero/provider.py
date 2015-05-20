from allauth.socialaccount import providers
from allauth.socialaccount.providers.base import (ProviderAccount,
                                                  AuthAction)
from allauth.socialaccount.providers.oauth.provider import OAuthProvider


class ZoteroAccount(ProviderAccount):
    pass

class ZoteroProvider(OAuthProvider):
    id = 'zotero'
    name = 'Zotero'
    package = 'users.providers.zotero'
    account_class = ZoteroAccount

    def get_default_scope(self):
        scope = []
        return scope

    def extract_uid(self, data):
        return data['userID']

    def extract_common_fields(self, data):
        return dict(username=data.get('username', ''))

providers.registry.register(ZoteroProvider)
