# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


from social.backends.oauth import BaseOAuth2
import logging
import requests
import orcid
from nameparser import HumanName
from .parsers import ParserOrcid
from django.db import transaction
from django.utils import timezone
from etalia.core.managers import PaperManager
from etalia.library.models import PaperUser
from .exceptions import OrcidInstanceError
from requests.exceptions import HTTPError

logger = logging.getLogger(__name__)


class OrcidOAuth2(BaseOAuth2):
    """Orcid OAuth authentication backend"""

    _type = 'ORC'
    parser = ParserOrcid()

    name = 'orcid'
    print_name = 'ORCiD'
    ID_KEY = 'orcid'
    AUTHORIZATION_URL = 'https://orcid.org/oauth/authorize'
    ACCESS_TOKEN_URL = 'https://orcid.org/oauth/token'
    ACCESS_TOKEN_METHOD = 'POST'
    DEFAULT_SCOPE = ['/authenticate']
    SCOPE_SEPARATOR = ','
    REDIRECT_STATE = False
    STATE_PARAMETER = False
    REFERENCE_MANAGER = False
    EXTRA_DATA = [
        ('name', 'name'),
        ('access_token', 'access_token'),
        ('refresh_token', 'refresh_token'),
        ('token_type', 'token_type'),
        ('expires_in', 'expires_in'),
        ('scope', 'scope'),
    ]

    def get_user_data(self, access_token, *args, **kwargs):
        """Loads user data from service"""
        return self.get_json(
            'https://pub.orcid.org/v1.2/{id}'.format(id='test'),
            headers={'Authorization': 'Bearer {0}'.format(access_token)}
        )

    def get_user_details(self, response):
        name = HumanName(response.get('name'))
        return {'first_name': name.first.strip(),
                'last_name': name.last.strip()}

    def get_user_public_profile(self, user):

        session = requests.Session()
        key, secret = self.get_key_and_secret()

        api = orcid.PublicAPI(key, secret, sandbox=False)
        # Get token
        res = session.post(self.ACCESS_TOKEN_URL,
                          data={'client_id': key,
                                'client_secret': secret,
                                'scope': '/read-public',
                                'grant_type': 'client_credentials'},
                          headers={'Accept': 'application/json'})
        res.raise_for_status()

        access_token = res.json().get('access_token')
        uid = user.social_auth.get(provider=self.name).uid
        res = session.get('https://pub.orcid.org/v1.2/{0}/orcid-profile/'.format(uid),
                           headers={
                               'Content-Type': 'application/orcid+json',
                               'Authorization': 'Bearer ' + access_token
                           })
        res.raise_for_status()

        return res.json()

    def get_public_papers(self, user):

        data = self.get_user_public_profile(user)
        logger.info(data)
        if data.get('orcid-profile').get('orcid-activities'):
            return data.get('orcid-profile')\
                .get('orcid-activities')\
                .get('orcid-works')\
                .get('orcid-work')
        else:
            return []

    def update_lib(self, user, full=False):

        try:
            papers = self.get_public_papers(user)
        except HTTPError as exc:
            raise OrcidInstanceError(user, exc.response)
        count = 0
        not_new_stack_count = 0

        for item in papers:

            try:
                entry = self.parser.parse(item)
            except Exception:
                logger.exception('Orcid parser failed')
                continue

            try:
                paper, journal = self.add_entry(entry)
            except Exception:
                logger.exception('Orcid adding paper failed')
                continue

            if paper:
                logger.info(
                    '+ Entry: {ids} from {user} / {backend}'.format(
                        ids=paper.print_ids,
                        user=user.email,
                        backend=self.name))
                orcid_id = item.get('source', {})\
                    .get('source-orcid', {})\
                    .get('path')
                info = {'authored': True,
                        'created': paper.date_pp or timezone.now().date()}
                new = self.associate_paper(user,
                                           paper,
                                           orcid_id,
                                           info)

                if new:
                    count += 1
                    not_new_stack_count = 0
                else:
                    not_new_stack_count += 1

            if not full:
                if not_new_stack_count > 50:
                    break  # exit when reaching a stack of 10 already uploaded references

        return count

    @staticmethod
    def associate_paper(user, paper, provider_id, info):
        """Update PaperUser and UserLibPaper table"""
        with transaction.atomic():
            pu, new = PaperUser.objects.get_or_create(user=user, paper=paper)
            pu.add(provider_id=provider_id, info=info, orcid=True)
        return new

    def add_entry(self, entry):
        """ConsolidateManager and add entry to DB"""

        # insert to DB or retrieve
        entry['is_trusted'] = False
        paper, journal = PaperManager(consolidate=True)\
            .get_or_create_from_entry(entry)

        return paper, journal