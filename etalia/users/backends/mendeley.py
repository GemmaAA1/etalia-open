# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import sys
import logging
from dateutil.parser import parse
from social.backends.oauth import BaseOAuth2
from social.backends.mendeley import MendeleyMixin
from mendeley import Mendeley
from mendeley.models.common import Person
from mendeley.auth import MendeleySession, \
    MendeleyAuthorizationCodeTokenRefresher
from mendeley.exception import MendeleyApiException
from .BaseMixin import BackendLibMixin
from .parsers import ParserMendeley
from time import time

logger = logging.getLogger(__name__)


class CustomMendeleyOAuth2(MendeleyMixin, BackendLibMixin, BaseOAuth2):

    # library specific
    _type = 'MEN'
    parser = ParserMendeley()

    # base
    name = 'mendeley'
    print_name = 'Mendeley'
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
        if 'employment' in response:
            out['tmp_affiliation'] = {}
            first_aff = response.get('employment')[0]
            first_aff_details = first_aff.get('institution_details', '')
            out['tmp_affiliation']['department'] = first_aff.get('department', '')
            out['tmp_affiliation']['institution'] = first_aff.get('institution', '')
            out['tmp_affiliation']['city'] = first_aff_details.get('city', '')
            out['tmp_affiliation']['state'] = first_aff_details.get('state', '')
            out['tmp_affiliation']['country'] = first_aff_details.get('country', '')
        try:
            out['title'] = response.get('title', '')
            out['position'] = response.get('academic_status', '')
            out['photo'] = response.get('photo', '')['standard']
        except KeyError:
            out['photo'] = ''
        return out

    def get_session(self, social, user, *args, **kwargs):

        client_id, client_secret = self.get_key_and_secret()
        tokens = {
            'access_token': social.access_token,
            'refresh_token': social.extra_data['refresh_token']
        }

        # start authorization flow
        mendeley = Mendeley(
            client_id,
            client_secret=client_secret,
            redirect_uri=self.redirect_uri
        )
        auth = mendeley.start_authorization_code_flow()
        refresher = MendeleyAuthorizationCodeTokenRefresher(auth)
        session = MendeleySession(
            mendeley,
            tokens,
            client=auth.client,
            refresher=MendeleyAuthorizationCodeTokenRefresher(auth)
        )

        # test token expiration
        expires_at = social.extra_data['expires_at']
        if (expires_at or 0) < time():   # renew
            refresher.refresh(session)

            # Store new tokens
            for key in session.token.keys():
                social.extra_data[key] = session.token[key]
            social.save(update_fields=('extra_data', ))

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

        # retrieve list of documents per page
        page = session.documents.list(
            page_size=self.CHUNK_SIZE,
            sort='created',
            order='desc',
            view='all')

        # loop through documents/page until end or find already registered doc
        while True:
            for item in page.items:

                try:
                    entry = self.parser.parse(item['data'])
                except Exception:
                    logger.exception(sys.exc_info())
                    continue

                try:
                    paper, journal = self.add_entry(entry)
                except Exception:
                    logger.error(sys.exc_info())
                    continue

                if paper:
                    logger.info(
                        '+ Entry: {ids} from {user} / {backend}'.format(
                            ids=paper.print_ids,
                            user=user.email,
                            backend=self.name))
                    new = self.associate_paper(user,
                                               paper,
                                               item.id,
                                               entry['user_info'])
                    if new:
                        count += 1
                        not_new_stack_count = 0
                    else:
                        not_new_stack_count += 1
                else:
                    logger.info(
                        '- Item: {type_} from {user} / {backend}'.format(
                            type_=item.type,
                            user=user.email,
                            backend=self.name))

            if not full:
                if not_new_stack_count > 50:
                    break  # exit when reaching 50 already uploaded references

            if page.next_page:
                page = page.next_page
            else:
                break

        return count

    @staticmethod
    def build_mendeley_identifiers(paper):
        identifiers = {}
        if paper.id_doi:
            identifiers['doi'] = paper.id_doi
        if paper.id_pii:
            identifiers['scopus'] = paper.id_pii
        if paper.id_pmi:
            identifiers['pmid'] = paper.id_pmi
        if paper.id_arx:
            identifiers['arxiv'] = paper.id_arx
        if paper.journal.id_issn or paper.journal.id_eissn:
            identifiers['issn'] = paper.journal.id_issn or paper.journal.id_eissn
        return identifiers

    def add_paper(self, session, paper):

        mend_doc_type = dict([(doctype[1], doctype[0])
                              for doctype in self.parser.MENDELEY_PT])
        if paper.type:
            type_ = mend_doc_type.get(paper.type, 'journal')  # 'PRE' (pre-print is unknow type for mendeley)
        else:
            type_ = 'journal'

        published_date = paper.date_ep or paper.date_pp or paper.date_fs

        authors = paper.authors.all()
        mend_authors = []
        for auth in authors:
            mend_authors.append(Person.create(auth.first_name, auth.last_name))
        try:
            ids = self.build_mendeley_identifiers(paper)
            resp = session.documents.create(
                paper.title,
                type_,
                identifiers=ids,
                websites=[paper.url],
                day=published_date.day,
                month=published_date.month,
                year=published_date.year,
                volume=paper.volume,
                issue=paper.issue,
                pages=paper.page,
                abstract=paper.abstract,
                source=paper.journal.title,
                authors=mend_authors,
            )
        except MendeleyApiException:
            return 1, None, None

        return None, resp.id, {'created': parse(str(resp.created) or 'Nothing'),
                               'last_modified': parse(str(resp.last_modified) or 'Nothing')}

    def trash_paper(self, session, paper_provider_id):
        try:
            doc = session.documents.get(id=paper_provider_id)
            if doc:
                resp = doc.move_to_trash()
                return 0
        except MendeleyApiException:
            return 1
