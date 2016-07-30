# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db.models import Q

from etalia.library.models import Paper, Journal, Author, AuthorPaper, CorpAuthor, \
    CorpAuthorPaper
from etalia.library.forms import PaperFormFillBlanks
from etalia.library.tasks import embed_paper

from ..constants import USERLIB_SYNCING, USERLIB_IDLE
from etalia.library.models import PaperUser


class BackendLibMixin(object):
    """Mixin for provider backend"""

    CHUNK_SIZE = 10
    _type = None
    parser = None

    def get_or_create_entry(self, entry):
        """Get or Create entry from/in Library

        Entry is a dictionary structure parsed by provider Parser.
        If corresponding paper is already in library (based on paper ids), the
        entry tried to consolidate the library

        Args:
            entry (dict): Entry to be added coming from consumer parser

        Returns:
            (Paper instance): Paper instance created
            (Journal instance): Corresponding Journal instance if found

        Raises:
        """
        try:
            # requirement to be a paper: have a title and an author and be a
            # supported paper type
            if entry['paper'].get('title', '') and entry['authors'] and \
                    entry['paper'].get('type', ''):
                item_paper = entry['paper']
                item_journal = entry['journal']
                item_paper['source'] = self._type
                # create/consolidate paper
                try:
                    paper = Paper.objects.get(Q(id_doi=item_paper['id_doi']) |
                                              Q(id_pmi=item_paper['id_pmi']) |
                                              Q(id_pii=item_paper['id_pii']) |
                                              Q(id_arx=item_paper['id_arx']) |
                                              Q(id_isbn=item_paper['id_isbn']) |
                                              Q(id_oth=item_paper['id_oth']))
                    if paper.is_trusted:
                        return paper, paper.journal

                except Paper.DoesNotExist:
                    paper = None

                form = PaperFormFillBlanks(item_paper, instance=paper)
                if form.is_valid():
                    paper = form.save()
                    # get journal
                    if self.is_journal_has_id(item_journal):
                        try:
                            # 1st and 2nd conditions are because most of providers
                            # do not distinguish e-issn and issn
                            journal = Journal.objects.get(
                                Q(id_issn=item_journal['id_issn']) |
                                Q(id_eissn=item_journal['id_issn']) |
                                Q(id_eissn=item_journal['id_eissn']) |
                                Q(id_arx=item_journal['id_arx']) |
                                Q(id_oth=item_journal['id_oth']))
                        except Journal.DoesNotExist:
                            journal = None
                    else:
                        try:
                            journal = Journal.objects.get(
                                title__iexact=item_journal['title'])
                        except Journal.DoesNotExist:
                            journal = None
                    paper.journal = journal
                    paper.is_trusted = False  # we do not trust provider source because they can be user made
                    paper.save()

                    # create/get authors
                    for pos, item_author in enumerate(entry['authors']):
                        author, _ = Author.objects.get_or_create(
                            first_name=item_author['first_name'],
                            last_name=item_author['last_name'])
                        AuthorPaper.objects.get_or_create(paper=paper,
                                                          author=author,
                                                          position=pos)
                    # create/get corp author
                    for pos, item_corp_author in enumerate(entry['corp_authors']):
                        corp_author, _ = CorpAuthor.objects.get_or_create(
                            name=item_corp_author['name']
                        )
                        CorpAuthorPaper.objects.get_or_create(
                            paper=paper,
                            corp_author=corp_author)
                    return paper, journal
            return None, None
        except Exception as e:
            return None, None

    def add_entry(self, entry):
        """Add or Retrieve entry and Update NLP Models"""

        paper, journal = self.get_or_create_entry(entry)

        if paper:
            embed_paper(paper.pk)

        return paper, journal

    @staticmethod
    def associate_paper(user, paper, provider_id, info):
        """Update PaperUser and UserLibPaper table"""
        pu, new = PaperUser.objects.get_or_create(user=user, paper=paper)
        pu.add(provider_id, info)
        # if new:
        #     pu.pin()
        return new

    def get_session(self, social, user, *args, **kwargs):
        raise NotImplementedError('Implement in subclass')

    def update_lib(self, user, session):
        # update db state
        user.lib.set_state(USERLIB_SYNCING)
        user.stats.log_lib_starts_sync(user)
        # really update
        count = self._update_lib(user, session)
        # retrieve first paper added
        user.lib.set_d_oldest()
        # update UserLib and Stats
        user.stats.log_lib_ends_sync(user, count)
        user.lib.set_state(USERLIB_IDLE)
        return count

    def _update_lib(self, session, user):
        raise NotImplementedError('Implement in subclass')

    def is_journal_has_id(self, item_journal):
        return any([True for key, val in item_journal.items()
                    if key.startswith('id_') and val])



