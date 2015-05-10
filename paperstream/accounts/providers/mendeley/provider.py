from allauth.socialaccount import providers
from allauth.socialaccount.providers.base import ProviderAccount
from allauth.socialaccount.providers.oauth2.provider import OAuth2Provider


class MendeleyAccount(ProviderAccount):

    def get_profile_url(self):
        return self.account.extra_data.get('link')

    def get_avatar_url(self):
        return self.account.extra_data.get('photos')[0]

    def to_str(self):
        dflt = super(MendeleyAccount, self).to_str()
        return self.account.extra_data.get('last_name', dflt)


class MendeleyProvider(OAuth2Provider):
    id = 'mendeley'
    name = 'Mendeley'
    package = 'paperstream.apps.accounts.providers.mendeley'
    account_class = MendeleyAccount

    def extract_extra_data(self, data):
        return data.get('data', {})

    def get_default_scope(self):
        return ['basic']

    def extract_uid(self, data):
        return str(data['id'])

    def extract_common_fields(self, data):
        return dict(first_name=data.get('first_name'),
                    last_name=data.get('last_name'),
                    email=data.get('email'))


providers.registry.register(MendeleyProvider)
