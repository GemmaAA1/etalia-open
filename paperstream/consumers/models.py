from django.db import models
from django.conf import settings
from django.utils import timezone

from Bio import Entrez
from Bio import Medline
from model_utils import Choices, fields

from library.models import Journal
from core.models import TimeStampedModel
from .constants import CONSUMER_TYPE


class Consumer(TimeStampedModel):

    # IDS
    type = models.CharField(max_length=4, choices=CONSUMER_TYPE)

    # name
    name = models.CharField(max_length=200)

    # number of papers retrieved per call
    ret_max = models.IntegerField(default=25)

    # number of past day to look through during initialization
    day0 = models.IntegerField(default=60)

    # Journal associated with consumer
    journals = models.ManyToManyField(Journal, through='ConsumerJournal')

    # TODO: Migration failed when class defined as abstract. Why is that??
    # class Meta:
    #     abstract = True


class PubmedConsumer(Consumer):

    # email
    email = settings.PUBMED_EMAIL

    @classmethod
    def create(cls):
        pubmed_consumer = cls(type='PUBM')
        return pubmed_consumer

    def consumes_journal(self, pk):
        """Consumes Pubmed API for journal with primary key pk
        """

        # Retrieve related ConsumerJournal
        try:
            cj = self.consumerjournal_set.get(journal=pk)
        except ConsumerJournal.DoesNotExist:
            raise ValueError('Journal (pk={0}) not linked to PubmedConsumer'
                             ''.format(pk))
        # Check status
        if not cj.status == 'idle':
            msg = 'ConsumerJournal status not idle ({0})'.format(cj.status)
            raise ValueError(msg)
        # Check if journal has ISSN
        if not cj.journal.id_issn and not cj.journal.id_eissn:
            msg = 'No id_issn or id_eissn for journal (pk={0})'.format(pk)
            raise ValueError(msg)
        else:
            issn = cj.journal.id_issn or cj.journal.id_eissn

        # Configure API
        # add email for contact
        Entrez.email = self.email

        # define query start date based on last scan or day0
        if cj.date_last_cons:
            start_date = cj.date_last_cons
        else:  # journal has never been scanned
            start_date = timezone.now() - timezone.timedelta(self.day0)

        # format date for API call
        start_date_q = start_date.strftime('%Y/%m/%d')

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
        # Parse items
        entries = Medline.parse(handle)



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
    STATUS = Choices('inactive', 'idle', 'in_queue', 'error')

    # journal
    journal = models.ForeignKey(Journal)

    # consumer
    consumer = models.ForeignKey(Consumer)

    # datetime of last consumption
    date_last_cons = models.DateTimeField(null=True, default=None)

    # Status monitor
    status = fields.StatusField()
    status_changed = fields.MonitorField(monitor='status')

    class Meta:
        unique_together = ('journal', 'consumer')