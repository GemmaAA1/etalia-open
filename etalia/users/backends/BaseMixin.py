# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import transaction
from etalia.core.managers import PaperManager
from etalia.library.models import PaperUser

from ..constants import USERLIB_SYNCING, USERLIB_IDLE


class BackendLibMixin(object):
    """Mixin for provider backend"""

    REFERENCE_MANAGER = True
    CHUNK_SIZE = 10
    _type = None
    parser = None

    def add_entry(self, entry):
        """ConsolidateManager and add entry to DB"""

        # insert to DB or retrieve
        entry['is_trusted'] = False
        paper, journal = PaperManager(consolidate=True)\
            .get_or_create_from_entry(entry)

        return paper, journal

    @staticmethod
    def associate_paper(user, paper, provider_id, info):
        """Update PaperUser and UserLibPaper table"""
        with transaction.atomic():
            pu, new = PaperUser.objects.get_or_create(user=user, paper=paper)
            pu.add(provider_id, info)
        return new

    def get_session(self, user):
        raise NotImplementedError('Implement in subclass')

    def update_lib(self, user, full=False):
        # update db state
        user.lib.set_state(USERLIB_SYNCING)
        user.stats.log_lib_starts_sync(user)
        # really update
        count = self._update_lib(user, full=full)
        # retrieve first paper added
        if user.lib.papers.count() > 0:
            user.lib.set_d_oldest()
        # update UserLib and Stats
        user.stats.log_lib_ends_sync(user, count)
        user.lib.set_state(USERLIB_IDLE)
        return count

    def _update_lib(self, user, full=None):
        raise NotImplementedError('Implement in subclass')

    def is_journal_has_id(self, item_journal):
        return any([True for key, val in item_journal.items()
                    if key.startswith('id_') and val])



