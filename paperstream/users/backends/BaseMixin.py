from django.db.models import Q
from django.utils import timezone
from library.models import Paper, Journal, Author, AuthorPaper, CorpAuthor, \
    CorpAuthorPaper
from library.forms import PaperFormFillBlanks

from ..models import UserLibPaper, UserLibJournal, UserStats

class BackendLibMixin(object):

    CHUNK_SIZE = 10
    _type = None
    parser = None

    def get_or_create_entry(self, entry):
        try:
            # minimum to be a paper: have a title and an author and be a support paper type
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
                        return paper

                except Paper.DoesNotExist:
                    paper = None

                form = PaperFormFillBlanks(item_paper, instance=paper)
                if form.is_valid():
                    paper = form.save()
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
                    paper.journal = journal
                    paper.is_trusted = False  # we trust consumer source only currently
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

    @staticmethod
    def associate_paper(paper, user, info):

        ulp, new = UserLibPaper.objects.get_or_create(userlib=user.lib,
                                                      paper=paper)
        ulp.date_created = info.get('created', None)
        ulp.last_date_modified = info.get('last_modified', None)
        ulp.authored = info.get('authored', None)
        ulp.starred = info.get('starred', None)
        ulp.scored = info.get('scored', 0.)
        ulp.save()
        return new

    @staticmethod
    def associate_journal(journal, user):
        ulj, new = UserLibJournal.objects.get_or_create(userlib=user.lib,
                                                        journal=journal)
        ulj.papers_in_journal += 1
        ulj.save()

    @staticmethod
    def create_lib_stats(user, count):
        UserStats.objects.create_lib_row(user, count)

    def get_session(self, social, user, *args, **kwargs):
        raise NotImplementedError('Implement in subclass')

    def update_lib(self, session, user):
        raise NotImplementedError('Implement in subclass')


