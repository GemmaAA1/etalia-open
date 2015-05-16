import logging

from django.db import models
from django.conf import settings
from django.utils import timezone
from django.db.models import Q

# from paperstream import capp

from Bio import Entrez
from Bio import Medline
from model_utils import Choices, fields

from library.models import Journal, AuthorPaper, Paper, Author, CorpAuthor, \
    CorpAuthorPaper
from core.models import TimeStampedModel
from .parsers import PubmedParser
from .constants import CONSUMER_TYPE

from library.forms import AuthorForm, PaperFormFillBlanks

logger = logging.getLogger(__name__)


class Consumer(TimeStampedModel):
    """Kind of abstract Consumer class (kind of only because I cannot make
    it abstract with the manytomany field
    """
    #TODO: How to make this class abstract? See comment below

    # IDS
    type = models.CharField(max_length=4, choices=CONSUMER_TYPE)

    # name
    name = models.CharField(max_length=200, unique=True)

    # number of papers retrieved per call
    ret_max = models.IntegerField(default=25)

    # number of past day to look through during initialization
    day0 = models.IntegerField(default=61)

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
        cj, new = ConsumerJournal.objects.get_or_create(journal=journal,
                                                        consumer=self)
        return cj

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
                errors.append(err)
                pass
        return errors

    def activate_all(self):
        """Deactivate all journals in ConsumerJournal table when feasible"""
        errors = []
        for journal in self.journals.all():
            try:
                self.activate_journal(journal)
            except ValueError as err:
                errors.append(err)
                pass
        return errors

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
        if cj.status == 'idle' and (journal.id_issn or journal.id_eissn):
            return True
        else:
            return False

    def consume_journal(self, journal):
        """ Consume Consumer API and update stats

        :param: journal (Journal): journal instance
        :return: list of entries
        """
        raise NotImplemented

    def populate_journal(self, journal):
        """Check journal validity, consume api, save stats, parse entries,
        save records to DB

        :param: journal (Journal): journal instance
        """
        if self.journal_is_valid(journal):
            # retrieve new entries from journal
            entries = self.consume_journal(journal)

            # save to database
            for entry in entries:
                item = self.parser.parse(entry)
                item_paper = item['paper']
                item_paper['source'] = self.type
                # create/consolidate paper
                try:
                    paper = Paper.objects.get(Q(id_doi=item_paper['id_doi']) |
                                              Q(id_pmi=item_paper['id_pmi']) |
                                              Q(id_pii=item_paper['id_pii']) |
                                              Q(id_arx=item_paper['id_arx']) |
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
                for pos, item_author in enumerate(item['authors']):
                    author, _ = Author.objects.get_or_create(
                        first_name=item_author['first_name'],
                        last_name=item_author['last_name'])
                    AuthorPaper.objects.get_or_create(paper=paper,
                                                      author=author,
                                                      position=pos)
                # create/get corp author
                for pos, item_corp_author in enumerate(item['corp_authors']):
                    corp_author, _ = CorpAuthor.objects.get_or_create(
                        name=item_corp_author['name']
                    )
                    CorpAuthorPaper.objects.get_or_create(
                        paper=paper,
                        corp_author=corp_author)


class ConsumerPubmed(Consumer):
    """Pubmed consumer subclass
    """

    def __init__(self, *args, **kwargs):
        super(ConsumerPubmed, self).__init__(*args, **kwargs)
        self.type = 'PUBM'

    parser = PubmedParser()

    # email
    email = settings.PUBMED_EMAIL

    # objects = ConsumerPubmedManager()

    def consume_journal(self, journal):
        """Consumes Pubmed API for journal

        :param journal (Journal):
        :return: (list) of entry
        """

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status = 'consuming'
        cj.save()

        # Configure API
        # add email for contact
        Entrez.email = self.email

        try:
            # define query start date based on last scan or day0
            if cj.last_date_cons:
                start_date = cj.last_date_cons
            else:  # journal has never been scanned
                start_date = timezone.now() - timezone.timedelta(self.day0)
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

            entries = []
            for record in records:
                entries.append(record)

            # close handle
            handle.close()

            # Update consumer_journal
            cj.update(True, len(id_list))

        except Exception as e:
            cj.update(False, 0)

            logger.exception(e)
            logger.warning()
            return list()

        return entries


class ElsevierConsumer(Consumer):

    # Type
    type = 'ELSV'

    # API key
    api_key = settings.ELSEVIER_API_KEY

    # URL
    URL_QUERY = 'http://api.elsevier.com/content/search/scidir?query='


class ArxivConsumer(Consumer):

    # Type
    type = 'ARXI'

    URL_QUERY = 'http://export.arxiv.org/api/query?search_query='


class ConsumerJournal(models.Model):

    #
    STATUS = Choices('inactive', 'idle', 'in_queue', 'consuming', 'error')

    # journal
    journal = models.ForeignKey(Journal)

    # consumer
    consumer = models.ForeignKey(Consumer)

    # last update
    # datetime of last consumption
    last_date_cons = models.DateTimeField(null=True, default=None)
    # number of papers
    last_number_papers = models.IntegerField(default=0)

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

    def update(self, success, n):
        """Update attributes and ConsumerJournalStats

        :param success (bool):
        :param n (int): number of papers fetched
        """
        if success:
            self.status = 'idle'
            self.save()
            self.last_date_cons = self.status_changed
            self.last_number_papers = n
            self.stats.create(
                date=self.last_date_cons,
                number_papers=self.last_number_papers,
                status='SUC')
            self.save()
        else:
            self.status = 'error'
            self.save()
            self.last_date_cons = self.status_changed
            self.last_number_papers = 0
            self.stats.create(date=self.last_date_cons, number_papers=0,
                              status='FAI')
            self.save()


class ConsumerJournalStat(models.Model):

    consumer_journal = models.ForeignKey(ConsumerJournal, related_name='stats')

    # date of consumption
    date = models.DateTimeField(null=False)

    # number of papers
    number_papers = models.IntegerField(default=0)

    # status
    status = models.CharField(max_length=3,
                              choices=(('SUC', 'Success'), (('FAI'), 'Failed')))