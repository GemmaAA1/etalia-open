# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from social.backends.oauth import BaseOAuth1
from pyzotero import zotero
import logging

from .BaseMixin import BackendLibMixin
from .parsers import ParserZotero

logger = logging.getLogger(__name__)


class CustomZoteroOAuth(BackendLibMixin, BaseOAuth1):
    """Zotero OAuth authorization mechanism"""

    # library specific
    _type = 'ZOT'
    parser = ParserZotero()

    # base
    name = 'custom-zotero'
    print_name = 'Zotero'
    REQUIRES_EMAIL_VALIDATION = True
    AUTHORIZATION_URL = 'https://www.zotero.org/oauth/authorize'
    REQUEST_TOKEN_URL = 'https://www.zotero.org/oauth/request'
    ACCESS_TOKEN_URL = 'https://www.zotero.org/oauth/access'

    def get_user_id(self, details, response):
        """Return user unique id provided by service"""
        return details['userID']

    def get_user_details(self, response):
        """Return user details from Zotero API account"""
        access_token = response.get('access_token', {})
        return {
            'username': access_token.get('username', ''),
            'userID': access_token.get('userID', '')
        }

    def get_session(self, social, user, *args, **kwargs):
        """Return OAuth session"""
        # Authenticate to zotero
        uid = social.tokens['userID']
        token = social.tokens['oauth_token']
        session = zotero.Zotero(uid, 'user', token)

        return session

    def update_lib(self, user, session):
        """Update User Lib

        Args:
            user (User): A User instance
            session: An OAuth session as provided by get_session()

        Returns:
            (int): Number of papers added
        """
        # Init
        new = True
        count = 0
        not_new_stack_count = 0

        # update db states
        user.stats.log_lib_starts_sync(user)
        user.lib.set_state('ING')

        items = session.top(limit=self.CHUNK_SIZE)

        while True:
            for item in items:
                try:
                    entry = self.parser.parse(item['data'])
                except Exception as e:
                    logger.exception('Zotero parser failed')
                    continue
                paper, journal = self.add_entry(entry)

                if paper:
                    logger.info(
                        '+ Entry: {ids} from {user} / {backend}'.format(
                            ids=paper.print_ids,
                            user=user.email,
                            backend=self.name))
                    new = self.associate_paper(paper, user, entry['user_info'])
                    if new:
                        count += 1
                        not_new_stack_count = 0
                    else:
                        not_new_stack_count += 1
                    if journal:
                        self.associate_journal(journal, user)
                else:
                    logger.info(
                        '- Item: {type_} from {user} / {backend}'.format(
                            type_=item['data']['itemType'],
                            user=user.email,
                            backend=self.name))
            if not_new_stack_count > 10:
                break  # exit when reaching a stack of 10 already uploaded references
            try:
                items = session.follow()
            except:
                break

        # update UserLib and Stats
        user.stats.log_lib_ends_sync(user, count)
        user.lib.set_state('IDL')
        return count
