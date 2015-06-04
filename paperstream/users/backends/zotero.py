from social.backends.oauth import BaseOAuth1
from pyzotero import zotero

from .BaseMixin import BackendLibMixin
from .parsers import ParserZotero


class CustomZoteroOAuth(BackendLibMixin, BaseOAuth1):

    """Zotero OAuth authorization mechanism"""

    # library specific
    _type = 'ZOT'
    parser = ParserZotero()

    # base
    name = 'custom-zotero'
    AUTHORIZATION_URL = 'https://www.zotero.org/oauth/authorize'
    REQUEST_TOKEN_URL = 'https://www.zotero.org/oauth/request'
    ACCESS_TOKEN_URL = 'https://www.zotero.org/oauth/access'

    def get_user_id(self, details, response):
        """
        Return user unique id provided by service. For Ubuntu One
        the nickname should be original.
        """
        return details['userID']

    def get_user_details(self, response):
        """Return user details from Zotero API account"""
        access_token = response.get('access_token', {})
        return {
            'username': access_token.get('username', ''),
            'userID': access_token.get('userID', '')
        }

    def get_session(self, social, user, *args, **kwargs):

        # Authenticate to zotero
        uid = social.tokens['userID']
        token = social.tokens['oauth_token']
        session = zotero.Zotero(uid, 'user', token)

        return session

    def update_lib(self, session, user, *args, **kwargs):

        # Init
        new = True
        count = 0
        user.lib.set_lib_syncing()

        items = session.top(limit=self.CHUNK_SIZE)

        while True:
            for item in items:
                entry = self.parser.parse(item['data'])
                paper, journal = self.add_entry(entry)
                if paper:
                    new = self.associate_paper(paper, user, entry['user_info']) and new
                    if new:
                        count += 1
                    else:
                        break
                    if journal:
                        self.associate_journal(journal, user)
            try:
                items = session.follow()
            except:
                break

        # update UserLib and Stats
        self.create_lib_stats(user, count)
        user.lib.set_lib_idle()
        return count