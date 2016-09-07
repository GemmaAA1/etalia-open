# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import sys
from dateutil.parser import parse
from social.backends.oauth import BaseOAuth1
from pyzotero import zotero
from pyzotero.zotero_errors import ResourceNotFound
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
    name = 'zotero'
    print_name = 'Zotero'
    REQUIRES_EMAIL_VALIDATION = True
    AUTHORIZATION_URL = 'https://www.zotero.org/oauth/authorize'
    REQUEST_TOKEN_URL = 'https://www.zotero.org/oauth/request'
    ACCESS_TOKEN_URL = 'https://www.zotero.org/oauth/access'
    GET = {'write_access': '1'}

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

    def _update_lib(self, user, session, full=False):
        """Update User Lib

        Args:
            user (User): A User instance
            session: An OAuth session as provided by get_session()
            full (bool): if True, run a full update on library. If false, only
                update new added added matches

        Returns:
            (int): Number of matches added
        """
        # Init
        new = True
        count = 0
        not_new_stack_count = 0

        items = session.top(limit=self.CHUNK_SIZE)

        while True:
            for item in items:

                try:
                    entry = self.parser.parse(item['data'])
                except Exception:
                    logger.exception('Zotero parser failed')
                    continue

                try:
                    paper, journal = self.add_entry(entry)
                except Exception:
                    logger.exception('Zotero adding paper failed')
                    continue

                if paper:
                    logger.info(
                        '+ Entry: {ids} from {user} / {backend}'.format(
                            ids=paper.print_ids,
                            user=user.email,
                            backend=self.name))
                    new = self.associate_paper(user,
                                               paper,
                                               item['key'],
                                               entry['user_info'])

                    if new:
                        count += 1
                        not_new_stack_count = 0
                    else:
                        not_new_stack_count += 1
                else:
                    logger.info(
                        '- Item: {type_} from {user} / {backend}'.format(
                            type_=item['data']['itemType'],
                            user=user.email,
                            backend=self.name))
            if not full:
                if not_new_stack_count > 50:
                    break  # exit when reaching a stack of 10 already uploaded references
            try:
                items = session.follow()
            except:
                break

        return count

    def add_paper(self, session, paper):

        zot_doc_type = dict([(doctype[1], doctype[0])
                             for doctype in self.parser.ZOTERO_PT])
        # if paper.type:
        #     type_ = zot_doc_type[paper.type]
        # else:
        type_ = 'journalArticle'

        template = session.item_template(type_)

        template['title'] = paper.title
        template['DOI'] = paper.id_doi
        template['url'] = paper.url
        template['date'] = (paper.date_ep or paper.date_pp or paper.date_fs).strftime('%m/%d/%Y')
        template['volume'] = paper.volume
        template['issue'] = paper.issue
        template['pages'] = paper.page
        template['publicationTitle'] = paper.journal.title
        if paper.journal.id_issn and paper.journal.id_eissn:
            template['ISSN'] = '{issn},{e_issn}'.format(issn=paper.journal.id_issn,
                                                        e_issn=paper.journal.id_eissn)
        elif paper.journal.id_issn:
            template['ISSN'] = paper.journal.id_issn
        elif paper.journal.id_eissn:
            template['ISSN'] = paper.journal.id_eissn
        else:
            template['ISSN'] = ''
        authors = paper.authors.all()
        template['creators'] = []
        for auth in authors:
            template['creators'].append({'creatorType': 'author',
                                         'firstName': auth.first_name,
                                         'lastName': auth.last_name})
        if paper.id_arx:
            template['libraryCatalog'] = 'arXiv.org'

        # push
        resp = session.create_items([template])
        if not resp['failed']:
            return None, \
                   resp['success']['0'], \
                   {'created': parse(str(resp['successful']['0']['data']['dateAdded'])),
                    'last_modified': parse(str(resp['successful']['0']['data']['dateModified']))
                    }
        else:
            logger.warning(resp['failed'])
            return None, None, None

    def trash_paper(self, session, paper_provider_id):
        """Trash item from zotero library"""
        try:  # retrieve item by id
            item = session.item(paper_provider_id)
            if session.delete_item(item):
                return 0
        except ResourceNotFound:
            return 1



