import re
import logging
import requests
import feedparser
import time
from dateutil import parser
from django.db import models
from django.conf import settings
from django.utils import timezone
from django.db.models import Q, F

from config.celery import celery_app as app
from celery.contrib.methods import task_method


from Bio import Entrez
from Bio import Medline
from model_utils import Choices, fields

from core.utils import get_env_variable
from library.models import Journal, AuthorPaper, Paper, Author, CorpAuthor, \
    CorpAuthorPaper
from core.models import TimeStampedModel
from .parsers import ParserPubmed, ParserArxiv, ParserElsevier
from .constants import CONSUMER_TYPE
from nlp.tasks import all_embeddings_and_neighbors

from library.forms import AuthorForm, PaperFormFillBlanks

logger = logging.getLogger(__name__)

class Consumer(TimeStampedModel):
    """Kind of abstract Consumer class (kind of only because I cannot make
    it abstract with the manytomany field
    """
    #TODO: How to make this class abstract? See comment below

    # IDS
    type = models.CharField(max_length=3, choices=CONSUMER_TYPE)

    # name
    name = models.CharField(max_length=200, unique=True)

    # number of papers retrieved per call
    ret_max = models.IntegerField(default=25)

    # number of past day to look through during initialization
    day0 = models.IntegerField(default=settings.CONS_INIT_PAST)

    # Journal associated with consumer
    journals = models.ManyToManyField(Journal, through='ConsumerJournal')

    parser = None

    # Class canNot be abstract because through need to be child dependent.
    # See Ticket #11760 (https://code.djangoproject.com/ticket/11760#no1)
    # class Meta:
    #     abstract = True

    def add_journal(self, journal):
        """Link journal to consumer
        :param: journal (Journal)
        :return: (ConsumerJournal)
        """
        cj, new = ConsumerJournal.objects.get_or_create(
            journal=journal,
            consumer=self,
            last_date_cons=timezone.now() - timezone.timedelta(days=self.day0))
        return cj

    def add_journals(self, journals):
        for journal in journals:
            self.add_journal(journal)

    def activate_journal(self, journal):
        """Activate row in ConsumerJournal table is feasible
        :param: journal (Journal): journal instance
        """
        try:
            self.consumerjournal_set.get(journal=journal).activate()
        except Journal.DoesNotExist:
            msg = 'Journal not associated with consumer ' \
                  '(title: {0})'.format(journal.title)
            raise ValueError(msg)

    def deactivate_journal(self, journal):
        """Deactivate row in ConsumerJournal table is feasible
        :param: journal (Journal): journal instance
        """
        try:
            self.consumerjournal_set.get(journal=journal).deactivate()
        except Journal.DoesNotExist:
            msg = 'Journal not associated with consumer ' \
                  '(title: {0})'.format(journal.title)
            raise ValueError(msg)

    def deactivate_all(self):
        """Deactivate all journals in ConsumerJournal table when feasible"""
        errors = []
        for journal in self.journals.all():
            try:
                self.deactivate_journal(journal)
            except ValueError or AssertionError as err:
                logger.warning(err)
                pass

    def activate_all(self):
        """Deactivate all journals in ConsumerJournal table when feasible"""
        errors = []
        for journal in self.journals.all():
            try:
                self.activate_journal(journal)
            except ValueError as err:
                logger.warning(err)
                pass

    def reset_last_date_cons(self, journal):
        q = self.consumerjournal_set.filter(journal=journal)
        q.update(last_date_cons=
                 timezone.now() - timezone.timedelta(days=self.day0))

    def reset_all_last_date_cons(self):
        ConsumerJournal.objects.filter(consumer=self).update(
            last_date_cons=timezone.now() - timezone.timedelta(days=self.day0))

    def journal_is_valid(self, journal):
        """ Check if journal is valid
        :param journal (Journal): journal instance
        :return: (Bool)
        """
        # exist in ConsumerJournal
        try:
            cj = self.consumerjournal_set.get(journal=journal)
        except ConsumerJournal.DoesNotExist:
            return False
        # is 'idle' and journal has id_eissn or id_issn
        if cj.status == 'idle':
            return True
        else:
            return False

    def consume_journal(self, journal):
        """ Consume Consumer API and update stats

        :param: journal (Journal): journal instance
        :return: list of entries
        """
        raise NotImplemented

    def add_entry(self, entry, journal):
        try:
            # minimum to be a paper: have a title and an author
            if entry['paper'].get('title', '') and entry['authors']:
                item_paper = entry['paper']
                item_paper['source'] = self.type
                # create/consolidate paper
                try:
                    paper = Paper.objects.get(Q(id_doi=item_paper['id_doi']) |
                                              Q(id_pmi=item_paper['id_pmi']) |
                                              Q(id_pii=item_paper['id_pii']) |
                                              Q(id_arx=item_paper['id_arx']) |
                                              Q(id_isbn=item_paper['id_isbn']) |
                                              Q(id_oth=item_paper['id_oth']))
                except Paper.DoesNotExist:
                    paper = None
                form = PaperFormFillBlanks(item_paper, instance=paper)
                if form.is_valid():
                    paper = form.save()
                    paper.journal = journal
                    paper.is_trusted = True  # we trust consumer source
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
                    return paper
                else:
                    return None
            return None
        # TODO: specify exception
        except Exception as e:
            return None

    def get_start_date(self, cj):
        if cj.last_date_cons:
            # for safety we look one day back
            start_date = cj.last_date_cons - timezone.timedelta(days=1)
        else:  # journal has never been scanned
            start_date = timezone.now() - timezone.timedelta(self.day0)
        return start_date

    @app.task(filter=task_method)
    def populate_journal(self, journal_pk):
        """Check journal validity, consume api, save stats, parse entries,
        save records to DB

        :param: journal_pk (Journal): pk of journal instance
        """

        journal = Journal.objects.get(pk=journal_pk)

        paper_added = 0

        if self.journal_is_valid(journal):
            # retrieve new entries from journal
            entries, success = self.consume_journal(journal)

            # save to database
            for entry in entries:
                item = self.parser.parse(entry)
                paper = self.add_entry(item, journal)
                if paper:
                    paper_added += 1

                    # Embed paper and get closest neighbors
                    all_embeddings_and_neighbors(paper.pk)

            # Update consumer_journal
            cj = self.consumerjournal_set.get(journal=journal)
            cj.update_stats(success, len(entries), paper_added)

            logger.info('populating {0}: ok'.format(journal.title))
        return paper_added

    def run_once_per_period(self):
        """
        This works as follows: run_once_per_period is designed to run
        periodically by a scheduler. To avoid consuming a bench of not
        very active journals (journals that published only every 3 months or
        so and there are a lot), we use a dynamic counter (<base_coundown_period>
        in [settings.CONS_MIN_DELAY, settings.CONS_MAX_DELAY]) that is updated
        at evry new consumption (increase if no paper fetched, decrease if
        paper fetched). After each consumption <coundown_period> is reset to
        <base_countdown_period>. Each new call to run_once_per_period decreases
        <coundown_period> by 1 and queue the ConsumerJournal instance for
        consumption is <countdown_period> < 1.

        :return:
        """

        logger.info('starting {0}:{1} daily consumption'.format(self.type,
                                                                self.name))

        # Get journal active
        consumer_journal_active = self.consumerjournal_set.filter(
            Q(status='idle') | Q(status='consuming') | Q(status='in_queue'))

        # Decrease day counter by 1
        consumer_journal_active.filter(coundown_period__lt=0).update(
            coundown_period=F('coundown_period')-1)

        # Get journal active and which counter at < 1
        consumerjournals_go_to_queue = \
            consumer_journal_active.filter(
                coundown_period__lt=1).select_related('journal')

        # queue journal for consumption
        for consumerjournal in consumerjournals_go_to_queue:
            self.populate_journal.apply_async(args=[consumerjournal.journal.pk, ],
                                              countdown=1,
                                              serializer='pickle')


class ConsumerPubmed(Consumer):
    """Pubmed consumer subclass
    """

    def __init__(self, *args, **kwargs):
        super(ConsumerPubmed, self).__init__(*args, **kwargs)
        self.type = 'PUB'

    parser = ParserPubmed()

    # email
    email = get_env_variable('PUBMED_EMAIL')

    def journal_is_valid(self, journal):
        if super(ConsumerPubmed, self).journal_is_valid(journal):
            if journal.id_issn or journal.id_eissn:
                return True

        return False

    def consume_journal(self, journal):
        """Consumes Pubmed API for journal

        :param journal (Journal):
        :return: (list) of entry
        """

        entries = []
        ok = True

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status = 'consuming'
        cj.save()

        # Configure API
        # add email for contact
        Entrez.email = self.email

        try:
            start_date = self.get_start_date(cj)
            # format date for API call
            start_date_q = start_date.strftime('%Y/%m/%d')

            issn = journal.id_issn or journal.id_eissn

            # build query (time range + journal issn)
            query = '(("{start_date}"[Date - Entrez] : "3000"[Date - Entrez])) AND ' \
                    '"{title}"[Journal]'.format(start_date=start_date_q,
                                                title=issn)

            # Call API
            # search items corresponding to query
            id_list = []
            ret_start = 0
            while True:
                # query
                handle = Entrez.esearch(db="pubmed",
                                        term=query,
                                        retmax=self.ret_max,
                                        retstart=ret_start)
                # read results
                record = Entrez.read(handle)
                # close handle
                handle.close()
                # get id list
                id_list = id_list + list(record["IdList"])
                if len(record["IdList"]) < int(self.ret_max):
                    break

                # update ret_start
                ret_start += self.ret_max
            # fetch items
            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline",
                                   retmode="text")

            records = Medline.parse(handle)

            for record in records:
                entries.append(record)

            # close handle
            handle.close()

            logger.debug('consuming {0}: OK'.format(journal.title))
        except Exception as e:
            ok = False
            entries = []
            logger.warning('consuming {0}: FAILED'.format(journal.title))

        return entries, ok


class ConsumerElsevier(Consumer):

    parser = ParserElsevier()

    # API key
    API_KEY = get_env_variable('ELSEVIER_API_KEY')

    # URL
    URL_QUERY = 'http://api.elsevier.com/content/search/index:SCIDIR?query='

    def __init__(self, *args, **kwargs):
        super(ConsumerElsevier, self).__init__(*args, **kwargs)
        self.type = 'ELS'

    def journal_is_valid(self, journal):
        if super(ConsumerElsevier, self).journal_is_valid(journal):
            if journal.id_issn or journal.id_eissn:
                return True

        return False

    def consume_journal(self, journal):
        # TODO: Volume is not fetched (
        """Consumes Elsevier API for journal

        ELSEVIER is a pain. I cannot make their API sorts the results by date
        accordingly to what is ahead of print and inprint. A possible reason
        is that prism:coverDisplayDate mixes up date of e-print and date of p-print.
        WARNING: Elsevier limits to 25 items / request.

        Approach followed here:
            - request results sorted by covDate
            - sort results by DOI (doi have timestamp embeded in them)
            - check if date for last retrieve results in prior last_cons_date
            - if not, keep going and increment starting point

        Other issues:
            Using SCOPUS doenot help because doesnot list aheadofprint papers.
            Using SCIDIR doesnot have the volume field (SCOPUS has it. wtf)

        :param journal (Journal):
        :return: (list) of entry
        """

        entries = []
        ok = True

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status = 'consuming'
        cj.save()

        try:
            headers = {'X-ELS-APIKey': self.API_KEY}
            query = '{url}ISSN({issn})&sort=-coverDate'.format(
                url=self.URL_QUERY,
                issn=journal.id_issn or journal.id_eissn,
            )

            count = 0

            start_date = self.get_start_date(cj)

            while True:
                data = {'count': str(self.ret_max),
                        'field': u'doi,coverDate,coverDisplayDate,url,identifier,title,'
                                 'publicationName,issueIdentifier,coverDisplayName,'
                                 'authors,creator,description,startingPage,issn,'
                                 'endingPage',
                        'start': str(count),
                        }

                resp = requests.post(query, data=data, headers=headers)

                entries += resp.json()['search-results']['entry']
                count += self.ret_max

                # sort entries by doi (doi are really time stamped)
                entries = sorted(entries, key=lambda x: x['prism:doi'][-11:],
                                 reverse=True)

                if entries:
                    current_start_date = entries[-1]['prism:coverDisplayDate']
                    # strip 'Available online' tag if in coverDisplayDate
                    current_start_date = re.sub(r'Available online', '',
                                                current_start_date).strip()
                    current_start_date = parser.parse(current_start_date)
                    current_start_date = current_start_date.replace(
                        tzinfo=timezone.pytz.timezone('UTC'))
                    if current_start_date < start_date:
                        break
                    else:
                        count += self.ret_max
                else:
                    break

            logger.debug('consuming {0}: OK'.format(journal.title))
        except Exception as e:
            ok = False
            entries - []
            logger.warning('consuming {0}: FAILED'.format(journal.title))

        return entries, ok


class ConsumerArxiv(Consumer):

    def __init__(self, *args, **kwargs):
        super(ConsumerArxiv, self).__init__(*args, **kwargs)
        self.type = 'ARX'

    parser = ParserArxiv()

    URL_QUERY = 'http://export.arxiv.org/api/query?search_query='

    def journal_is_valid(self, journal):
        if super(ConsumerArxiv, self).journal_is_valid(journal):
            if journal.id_arx:
                return True

        return False

    def consume_journal(self, journal):

        entries = []
        ok = True

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status = 'consuming'
        cj.save()
        try:
            # get category
            cat = re.sub(r'arxiv\.', '', cj.journal.id_arx)

            start_date = self.get_start_date(cj)

            # format date for API call
            start_date_q = '{year}{month:02d}{day:02d}'.format(
                year=start_date.year,
                month=start_date.month,
                day=start_date.day)

            end_date_q = '{year}{month:02d}{day:02d}'.format(
                year=timezone.now().year,
                month=timezone.now().month,
                day=timezone.now().day)

            count = 0
            total_entries = 1
            while count == 0 or count < total_entries:
                query = '{url}cat:{cat}+AND+submittedDate:[{start}+TO+{end}]' \
                        '&start={count}&max_results={ret_max}'.format(
                            url=self.URL_QUERY,
                            cat=cat,
                            start=start_date_q,
                            end=end_date_q,
                            count=count,
                            ret_max=self.ret_max)
                resp = requests.get(query)
                time.sleep(1)
                data = feedparser.parse(resp.text)
                total_entries = int(data['feed']['opensearch_totalresults'])
                if len(data.entries) < 25 and \
                                count < (total_entries - self.ret_max):
                    # retry once
                    resp = requests.get(query)
                    time.sleep(1)  # for politeness
                    data = feedparser.parse(resp.text)
                count += self.ret_max
                entries += data.entries

            logger.debug('consuming {0}: OK'.format(journal.title))
        except Exception as e:
            ok = False
            entries = []
            logger.warning('consuming {0}: FAILED'.format(journal.title))

        return entries, ok


class ConsumerJournal(models.Model):

    #
    STATUS = Choices('inactive', 'idle', 'in_queue', 'consuming', 'error')

    # journal
    journal = models.ForeignKey(Journal)

    # consumer
    consumer = models.ForeignKey(Consumer)

    # last update
    # datetime of last consumption
    last_date_cons = models.DateTimeField(null=True, blank=True)
    # number of papers recorded in db
    last_number_papers_recorded = models.IntegerField(default=0)
    # number of papers fetched from provider
    last_number_papers_fetched = models.IntegerField(default=0)

    base_coundown_period = models.IntegerField(default=1)
    coundown_period = models.IntegerField(default=0)

    # Status monitor
    status = fields.StatusField(default='inactive')
    status_changed = fields.MonitorField(monitor='status')

    class Meta:
        unique_together = ('journal', 'consumer')

    def activate(self):
        if self.status == 'inactive':
            self.status = 'idle'
            self.save()
        elif self.status in ['idle', 'in_queue', 'consuming']:
            msg = 'journal already active: current status is {0}'.format(self.status)
            logger.info(msg)
        else:
            msg = 'cannot activate, journal status is {0}'.format(self.status)
            logger.info(msg)
            raise ValueError(msg)

    def deactivate(self):
        if self.status in ['idle', 'error']:
            self.status = 'inactive'
            self.save()
        elif self.status in ['in_queue', 'consuming']:
            msg = 'cannot deactivate, journal is busy: current status is {0}'.format(self.status)
            logger.info(msg)
            raise ValueError(msg)
        else:
            msg = 'journal already active: current status is {0}'.format(self.status)
            logger.info(msg)
            raise ValueError(msg)

    def update_stats(self, success, n_fet, n_rec):
        """Update attributes and ConsumerJournalStats

        :param success (bool):
        :param n_ret (int): number of papers fetched from API
        :param n_rec (int): number of papers recorded in db
        """
        if success:
            self.status = 'idle'
            self.last_date_cons = timezone.now()
            self.last_number_papers_fetched = n_fet
            self.last_number_papers_recorded = n_rec
            self.stats.create(
                  number_papers_fetched=self.last_number_papers_fetched,
                  number_papers_recorded=self.last_number_papers_recorded,
                  status='SUC')
            # update base countdown
            if n_fet < 0:
                if self.base_coundown_period < settings.CONS_MAX_DELAY:
                    self.base_coundown_period += 1
                else:
                    self.base_coundown_period = settings.CONS_MAX_DELAY
            elif n_fet > 0:
                if self.base_coundown_period > settings.CONS_MIN_DELAY:
                    self.base_coundown_period -= 1
                else:
                    self.base_coundown_period = settings.CONS_MIN_DELAY
            # reinit counter
            self.coundown_period = self.base_coundown_period
            self.save()
        else:
            self.status = 'error'
            self.stats.create(number_papers_fetched=0,
                              number_papers_recorded=0,
                              status='FAI')
            self.save()

    class Meta:
        ordering = ['last_date_cons']


class ConsumerJournalStat(models.Model):

    consumer_journal = models.ForeignKey(ConsumerJournal, related_name='stats')

    # date of consumption
    datetime = models.DateTimeField(null=False, auto_now_add=True)

    # number of papers fetched
    number_papers_fetched = models.IntegerField()

    # number of papers recorded
    number_papers_recorded = models.IntegerField()

    # status
    status = models.CharField(max_length=3,
                              choices=(('SUC', 'Success'),
                                       (('FAI'), 'Failed')))
