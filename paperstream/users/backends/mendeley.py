# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging
from social.backends.oauth import BaseOAuth2
from social.backends.mendeley import MendeleyMixin
from mendeley import Mendeley
from mendeley.models.common import Person
from mendeley.auth import MendeleySession, \
    MendeleyAuthorizationCodeTokenRefresher
from mendeley.exception import MendeleyApiException

from ..constants import MENDELEY_PT
from .BaseMixin import BackendLibMixin
from .parsers import ParserMendeley

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

    def update_lib(self, user, session):

        # Init
        new = True
        count = 0
        not_new_stack_count = 0

        # update db state
        user.stats.log_lib_starts_sync(user)
        user.lib.set_state('ING')

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
                    entry = self.parser.parse(item)
                except Exception as e:
                    logger.exception('Mendeley parser failed')
                    continue
                paper, journal = self.add_entry(entry)

                if paper:
                    logger.info(
                        '+ Entry: {ids} from {user} / {backend}'.format(
                            ids=paper.print_ids,
                            user=user.email,
                            backend=self.name))
                    new = self.associate_paper(paper, user, entry['user_info'],
                                               item.id) and new
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
                            type_=item.type,
                            user=user.email,
                            backend=self.name))

            if not_new_stack_count > 10:
                break  # exit when reaching 10 already uploaded references

            if page.next_page:
                page = page.next_page
            else:
                break

        # update UserLib and Stats
        user.stats.log_lib_ends_sync(user, count)
        user.lib.set_state('IDL')
        return count

    @staticmethod
    def add_paper(session, paper):

        mend_doc_type = dict([(doctype[1], doctype[0]) for doctype in MENDELEY_PT])
        if paper.type:
            type_ = mend_doc_type[paper.type]
        else:
            type_ = 'journal'

        published_date = paper.date_ep or paper.date_pp or paper.date_fs

        authors = paper.authors.all()
        mend_authors = []
        for auth in authors:
            mend_authors.append(Person.create(auth.first_name, auth.last_name))
        try:
            ids = paper.build_mendeley_identifiers()
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
                authors=mend_authors
            )
        except MendeleyApiException:
            return 1

        return None, resp.id

    @staticmethod
    def trash_paper(session, ulp):
        try:
            doc = session.documents.get(id=ulp.paper_provider_id)
            if doc:
                resp = doc.move_to_trash()
                return 0
        except MendeleyApiException:
            return 1






