# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import logging
import time

import requests
import feedparser
from dateutil import parser
from django.db import models
from django.conf import settings
from django.utils import timezone
from django.db.models import Q, F
from Bio import Entrez
from Bio import Medline
from model_utils import Choices, fields

from etalia.library.models import Journal, AuthorPaper, Paper, Author, \
    CorpAuthor, CorpAuthorPaper
from etalia.library.forms import PaperFormFillBlanks
from etalia.core.models import TimeStampedModel
from etalia.library.tasks import embed_paper

from .parsers import PubmedParser, ArxivParser, ElsevierParser
from .constants import CONSUMER_TYPE

logger = logging.getLogger(__name__)


class Consumer(TimeStampedModel):
    """ Abstract Consumer Table

    Consumers consumes matches from the web and populate the library
    """

    # IDS
    type = models.CharField(max_length=3, choices=CONSUMER_TYPE)

    # name
    name = models.CharField(max_length=200, unique=True)

    # number of matches retrieved per call
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

    def __str__(self):
        return self.name

    def add_journal(self, journal):
        """Add journal to consumer"""
        cj, new = ConsumerJournal.objects.get_or_create(journal=journal,
                                                        consumer=self)
        if new:
            cj.last_date_cons = timezone.now() - timezone.timedelta(days=self.day0)
            cj.save()
        return cj

    def add_journals(self, journals):
        """Add journals to consumer"""
        for journal in journals:
            self.add_journal(journal)

    def activate_journal(self, journal):
        """Activate journal if possible"""
        try:
            self.consumerjournal_set.get(journal=journal).activate()
        except Journal.DoesNotExist:
            msg = 'Journal not associated with consumer ' \
                  '(title: {0})'.format(journal.title)
            raise ValueError(msg)

    def deactivate_journal(self, journal):
        """Deactivate journal"""
        try:
            self.consumerjournal_set.get(journal=journal).deactivate()
        except Journal.DoesNotExist:
            msg = 'Journal not associated with consumer ' \
                  '(title: {0})'.format(journal.title)
            raise ValueError(msg)

    def deactivate_all(self):
        """Deactivate all journals if idle, report others"""
        self.consumerjournal_set\
            .filter(status__in=['idle'])\
            .update(status='inactive')
        logger.info('({0}) {1}: Deactivate all journals'.format(self.pk,
                                                                self.name))
        cjs_failed = self.consumerjournal_set\
            .filter(status__in=['error', 'in_queue', 'consuming'])\
            .values('journal__pk', 'journal__short_title', 'status')

        for cj in cjs_failed:
            err = '({0}) {1}: Deactivate ({2}) {3} FAILED, status is {4}'\
                .format(self.pk,
                        self.name,
                        cj['journal__pk'],
                        cj['journal__short_title'],
                        cj['status'])
            logger.info(err)

        return cjs_failed

    def activate_all(self):
        """Activate all journals if inactive, report errors"""
        self.consumerjournal_set\
            .filter(status__in=['inactive'])\
            .update(status='idle')
        logger.info('({0}) {1}: Activate all journals'.format(self.pk,
                                                              self.name))
        cjs_failed = self.consumerjournal_set\
            .filter(status__in=['error'])\
            .values('journal__pk', 'journal__short_title', 'status')

        for cj in cjs_failed:
            err = '({0}) {1}: Activate ({2}){3} FAILED, status is {4}'\
                .format(self.pk,
                        self.name,
                        cj['journal__pk'],
                        cj['journal__short_title'],
                        cj['status'])
            logger.info(err)

        return cjs_failed

    def reset_last_date_cons(self, journal):
        """Reset last date of consumption to default"""
        q = self.consumerjournal_set.filter(journal=journal)
        q.update(last_date_cons=timezone.now() -
                                timezone.timedelta(days=self.day0))

    def reset_all_last_date_cons(self):
        """Reset all last date of consumption to default"""
        ConsumerJournal.objects.filter(consumer=self).update(
            last_date_cons=timezone.now() - timezone.timedelta(days=self.day0))

    def journal_is_valid(self, journal):
        """ Check if journal is valid

        Args:
            journal: Journal instance

        Returns:
            (Bool): True if valid
        """
        # exist in ConsumerJournal
        try:
            cj = self.consumerjournal_set.get(journal=journal)
        except ConsumerJournal.DoesNotExist:
            return False
        # is 'idle' and journal has id_eissn or id_issn
        if cj.status == 'idle' or cj.status == 'retry':
            return True
        else:
            return False

    def consume_journal(self, journal):
        """ Consume Consumer API and update stats

        Args:
            journal: Journal instance

        Returns:
            (list): List of consumed entries
        """
        raise NotImplemented

    def add_entry(self, entry, journal):
        """Add entry to database library

        Entry is a dictionary structure parsed by consumer Parser.
        If corresponding paper is already in library (based on paper ids), the
        entry tried to consolidate the library paper

        Args:
            entry (dict): Entry to be added coming from consumer parser
            journal: Journal instance

        Returns:
            (Paper instance): Paper instance created
        """
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
                    # we trust consumer as source
                    paper.is_trusted = True
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
        except Exception:
            raise

    def get_start_date(self, cj):
        if cj.last_date_cons:
            # for safety we look one day back
            start_date = cj.last_date_cons - timezone.timedelta(days=1)
        else:  # journal has never been scanned
            start_date = timezone.now() - timezone.timedelta(self.day0)
        return start_date

    def populate_journal(self, journal_pk):
        """Consume data from journal

        Check journal validity, consume api, save stats, parse entries,
        and save records to database.

        Args:
            journal_pk (Journal): pk of journal instance
        """

        journal = Journal.objects.get(pk=journal_pk)

        paper_added = 0

        if self.journal_is_valid(journal):
            try:

                logger.info('Consuming {type}/{consumer}/{title} - starting...'
                    .format(type=self.type,
                            consumer=self.name,
                            title=journal.title))

                # retrieve new entries from journal
                entries = self.consume_journal(journal)

                # save to database
                for entry in entries:
                    item = self.parser.parse(entry)
                    paper = self.add_entry(item, journal)
                    if paper:
                        paper_added += 1

                        # Embed paper
                        embed_paper(paper.pk)

                # Update consumer_journal
                cj = self.consumerjournal_set.get(journal=journal)
                cj.update_stats('pass', len(entries), paper_added)

                logger.info(
                    'Consuming {type}/{consumer}/{title} - ({count}) done'.format(
                        type=self.type,
                        consumer=self.name,
                        title=journal.title,
                        count=paper_added)
                )
            except Exception:

                raise

    def run_once_per_period(self):
        """Run Consumer: Consumes active Journal associated with consumer.

        It works as follows: run_once_per_period beats periodically (every
        <period>). To avoid consuming a bench of not very active journals, time
        between 2 consumptions is dynamically adapted per journal.

        <base_countdown_period> defined the period between two consumptions
        (in <period> unit).It take values in range [settings.CONS_MIN_DELAY,
        settings.CONS_MAX_DELAY].
        <countdown_period> is init to <base_countdown_period> after journal
        has been consumed. It is decreased by 1 at each run_once_per_period call
        journal is queued for consumption if <countdown_period> = 1

        After consumption, <base_counter_period> is increased by 1 if no paper
        was fetched, decreased by 1 if papers were fetched.
        """
        from .tasks import populate_journal

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
        # NB: concurrency of consumer queue is 4 and we gently want to respect
        # a 1/s request.
        for sec, consumerjournal in enumerate(consumerjournals_go_to_queue):
            populate_journal.apply_async(
                args=[self, consumerjournal.journal.pk, ],
                countdown=sec * 1.1)


class ConsumerPubmed(Consumer):
    """Pubmed Consumer"""

    def __init__(self, *args, **kwargs):
        super(ConsumerPubmed, self).__init__(*args, **kwargs)
        self.type = 'PUB'

    parser = PubmedParser()

    # email
    email = settings.CONSUMER_PUBMED_EMAIL

    def journal_is_valid(self, journal):
        if super(ConsumerPubmed, self).journal_is_valid(journal):
            if journal.id_issn or journal.id_eissn:
                return True

        return False

    def consume_journal(self, journal):
        """Consumes Pubmed API for journal"""

        entries = []

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status_to('consuming')

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
        except Exception:
            cj.status_to('error')
            raise
        return entries


class ConsumerElsevier(Consumer):
    """Pubmed Consumer"""

    parser = ElsevierParser()

    # API key
    API_KEY = settings.CONSUMER_ELSEVIER_API_KEY

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
        # TODO: Volume is not fetched
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
            Using SCOPUS doenot help because doesnot list aheadofprint matches.
            Using SCIDIR doesnot have the volume field (SCOPUS has it. wtf)

        :param journal (Journal):
        :return: (list) of entry
        """

        entries = []

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status_to('consuming')

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
                    current_start_date = entries[-1]['prism:coverDate'][0].get('$')
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
        except Exception:
            cj.status_to('error')
            raise
        return entries


class ConsumerArxiv(Consumer):
    """Arxiv Consumer"""
    def __init__(self, *args, **kwargs):
        super(ConsumerArxiv, self).__init__(*args, **kwargs)
        self.type = 'ARX'

    parser = ArxivParser()

    URL_QUERY = 'http://export.arxiv.org/api/query?search_query='

    def journal_is_valid(self, journal):
        if super(ConsumerArxiv, self).journal_is_valid(journal):
            if journal.id_arx:
                return True

        return False

    def consume_journal(self, journal):
        """Consume Arxiv 'journal'"""
        entries = []

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status_to('consuming')
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
                time.sleep(1)  # for politeness
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
        except Exception:
            cj.status_to('error')
            raise

        return entries


class ConsumerJournal(TimeStampedModel):
    """Table for Consumer-Journal relationship"""

    STATUS = Choices('inactive',
                     'idle',
                     'in_queue',
                     'consuming',
                     'error',
                     'retry')

    # journal
    journal = models.ForeignKey(Journal)

    # consumer
    consumer = models.ForeignKey(Consumer)

    # last update
    # datetime of last consumption
    last_date_cons = models.DateTimeField(null=True, blank=True, default=None)
    # number of matches recorded in db
    last_number_papers_recorded = models.IntegerField(default=0)
    # number of matches fetched from provider
    last_number_papers_fetched = models.IntegerField(default=0)

    base_coundown_period = models.IntegerField(default=1)
    coundown_period = models.IntegerField(default=0)

    # Status monitor
    status = fields.StatusField(default='inactive')
    status_changed = fields.MonitorField(monitor='status')

    class Meta:
        unique_together = ('journal', 'consumer')
        ordering = ['last_date_cons']

    def __str__(self):
        return '{0}@{1}'.format(self.journal.short_title or self.journal.title,
                                self.consumer.name)

    def status_to(self, status):
        self.status = status
        self.save(update_fields=['status', ])

    def reset(self):
        fields_to_reset = ['status',
                           'last_date_cons',
                           'last_number_papers_fetched',
                           'last_number_papers_recorded',
                           'base_coundown_period']
        for field in fields_to_reset:
            setattr(self, field, self._meta.get_field(field).default)
        self.save()
        self.stats.create(status='RES',
                          message='Reset ConsumerJournal')

    def activate(self):
        if self.status == 'inactive':
            self.status = 'idle'
            self.save()
            msg = '({0}) {1}: Activate {2}'.format(self.consumer.pk,
                                                   self.consumer.name,
                                                   self.journal)
            logger.info(msg)
            self.stats.create(status='ACT')
        elif self.status == 'error':
            msg = '({0}) {1}: Activate {2} FAILED, status is {3}'\
                .format(self.consumer.pk,
                        self.consumer.name,
                        self.journal,
                        self.status)
            self.stats.create(status='FAI', message=msg)
            raise ValueError(msg)

    def deactivate(self):
        if self.status in ['idle', 'error']:
            self.status = 'inactive'
            self.save()
            msg = '({0}) {1}: Deactivate {2}'.format(self.consumer.pk,
                                                     self.consumer.name,
                                                     self.journal)
            logger.info(msg)
        elif self.status in ['in_queue', 'consuming']:
            msg = '({0}) {1}: Deactivate {2} FAILED, status is {3}'\
                .format(self.consumer.pk,
                        self.consumer.name,
                        self.journal,
                        self.status)
            logger.info(msg)
            self.stats.create(status='FAI', message=msg)
            raise ValueError(msg)

    def update_stats(self, success, n_fet, n_rec):
        """Update attributes and ConsumerJournalStats

        Args:
            success (bool): True if comsumption was a success
            n_ret (int): number of matches fetched from API
            n_rec (int): number of matches recorded in db
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

    def print_stats(self):
        tmp = ['Date\tState\t# Fetched\t# Recorded\tMessages\n']
        for stat in self.stats.all():
            tmp.append('{date}\t{state}\t{fetch}\t{reco}\t{mess}\n'.format(
                date=stat.datetime,
                state=stat.status,
                fetch=stat.number_papers_fetched,
                reco=stat.number_papers_recorded,
                mess=stat.message,
            ))
        print(''.join(tmp))


class ConsumerJournalStat(TimeStampedModel):
    """Table for ConsumerJournal stats"""
    consumer_journal = models.ForeignKey(ConsumerJournal, related_name='stats')

    # date of consumption
    datetime = models.DateTimeField(null=False, auto_now_add=True)

    # state
    status = models.CharField(max_length=3,
                              choices=(('SUC', 'Success'),
                                       ('FAI', 'Failed'),
                                       ('RES', 'Reset'),
                                       ('ACT', 'Activate'),
                                       ('DEA', 'Deactivate')))

    # number of matches fetched
    number_papers_fetched = models.IntegerField(default=0)

    # number of matches recorded
    number_papers_recorded = models.IntegerField(default=0)

    # message
    message = models.CharField(max_length=512, default='')

    class Meta:
        ordering = ['datetime']
