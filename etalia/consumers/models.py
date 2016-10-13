# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import re
import logging
import time
import requests
from requests import adapters
import json
import feedparser
from dateutil.parser import parse
from dateutil import parser
from django.db import models
from django.conf import settings
from django.utils import timezone
from django.db.models import Q, F
from Bio import Entrez
from Bio import Medline
from model_utils import Choices, fields

from config.celery import celery_app as app
from celery.contrib.methods import task_method

from etalia.library.models import Journal, AuthorPaper, Paper, Author, \
    CorpAuthor, CorpAuthorPaper
from etalia.threads.models import Thread, PubPeer, PubPeerComment
from etalia.core.models import TimeStampedModel

from .parsers import PubmedPaperParser, ArxivPaperParser, ElsevierPaperParser, \
    PubPeerThreadParser, BiorxivPaperParser, SpringerPaperParser
from etalia.core.managers import PaperManager, PubPeerManager
from .constants import CONSUMER_TYPE
from .mixins import SimpleCrawlerListMixin

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
    day0 = models.IntegerField(default=settings.CONSUMER_INIT_PAST)

    # Journal associated with consumer
    journals = models.ManyToManyField(Journal, through='ConsumerJournal')

    TYPE = None
    parser = None

    # Class canNot be abstract because through need to be child dependent.
    # See Ticket #11760 (https://code.djangoproject.com/ticket/11760#no1)
    # class Meta:
    #     abstract = True

    def __init__(self, *args, **kwargs):
        super(Consumer, self).__init__(*args, **kwargs)
        if not self.type:
            self.type = self.TYPE
        # attach session
        self.session = requests.Session()
        adapter = adapters.HTTPAdapter(max_retries=3)
        self.session.mount('http://', adapter)

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

        Args:
            entry (dict): dictionary structure parsed by consumer PaperParser.
            journal: Journal instance

        Returns:
            (Paper instance): Paper instance created
        """
        try:
            entry['is_trusted'] = True
            paper, _ = PaperManager().get_or_create_from_entry(
                entry, fetch_journal=False)
            if paper:
                paper.journal = journal
                paper.save(update_fields=['journal'])
            return paper
        except Exception:
            raise

    def get_start_date(self, cj):
        if cj.last_date_cons:
            # for safety we look one day back
            start_date = cj.last_date_cons - timezone.timedelta(days=1)
        else:  # journal has never been scanned
            start_date = timezone.now() - timezone.timedelta(self.day0)
        return start_date.date()

    def populate_journal(self, journal_pk):
        """Consume data from journal

        Check journal validity, consume api, save stats, parse entries,
        and save records to database.

        Args:
            journal_pk (Journal): pk of journal instance
        """

        journal = Journal.objects.get(pk=journal_pk)

        paper_count = 0

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
                        paper_count += 1

                # Update consumer_journal
                cj = self.consumerjournal_set.get(journal=journal)
                cj.update_stats('pass', len(entries), paper_count)

                logger.info(
                    'Consuming {type}/{consumer}/{title} - ({count}) done'.format(
                        type=self.type,
                        consumer=self.name,
                        title=journal.title,
                        count=paper_count)
                )
            except Exception:
                raise

    def run_once_per_period(self):
        """Run Consumer: Consumes active Journal associated with consumer.

        It works as follows: run_once_per_period beats periodically (every
        <period>). To avoid consuming a bench of not very active journals, time
        between 2 consumptions is dynamically adapted per journal.

        <base_countdown_period> defined the period between two consumptions
        (in <period> unit).It take values in range [settings.CONSUMER_MIN_DELAY,
        settings.CONSUMER_MAX_DELAY].
        <countdown_period> is init to <base_countdown_period> after journal
        has been consumed. It is decreased by 1 at each run_once_per_period call
        journal is queued for consumption if <countdown_period> = 1

        After consumption, <base_counter_period> is increased by 1 if no paper
        was fetched, decreased by 1 if papers were fetched.
        """
        from .tasks import populate_journal as populate_journal_async

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
            populate_journal_async.apply_async(
                args=[self.id, consumerjournal.journal.pk, ],
                countdown=sec * 1.1)


class ConsumerPubmed(Consumer):
    """Pubmed Consumer"""

    TYPE = 'PUB'
    parser = PubmedPaperParser()
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
                handle = Entrez.esearch("pubmed",
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
            if id_list:
                handle = Entrez.efetch("pubmed", id=id_list, rettype="medline",
                                       retmode="text")

                records = Medline.parse(handle)

                for record in records:
                    entries.append(record)

                # close handle
                handle.close()
        except Exception:
            cj.status_to('error')
            raise
        cj.status_to('idle')
        return entries


class ConsumerElsevier(Consumer):
    """Pubmed Consumer"""

    TYPE = 'ELS'
    parser = ElsevierPaperParser()
    API_KEY = settings.CONSUMER_ELSEVIER_API_KEY
    URL_QUERY = 'http://api.elsevier.com/content/search/index:SCIDIR?query='

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

                resp = self.session.post(query, data=data, headers=headers)
                if 'search-results' in resp.json().keys():
                    entries += resp.json().get('search-results').get('entry')
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
                        if current_start_date.date() < start_date:
                            break
                        else:
                            count += self.ret_max
                    else:
                        break
                else:
                    break
        except Exception:
            cj.status_to('error')
            raise
        cj.status_to('idle')
        return entries


class ConsumerArxiv(Consumer):
    """Arxiv Consumer"""

    TYPE = 'ARX'
    parser = ArxivPaperParser()
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
                resp = self.session.get(query)
                time.sleep(1)  # for politeness
                data = feedparser.parse(resp.text)
                total_entries = int(data['feed']['opensearch_totalresults'])
                if len(data.entries) < 25 and \
                                count < (total_entries - self.ret_max):
                    # retry once
                    resp = self.session.get(query)
                    time.sleep(1)  # for politeness
                    data = feedparser.parse(resp.text)
                count += self.ret_max
                entries += data.entries
        except Exception:
            cj.status_to('error')
            raise
        cj.status_to('idle')
        return entries


class ConsumerSpringer(Consumer):

    TYPE = 'SPR'
    parser = SpringerPaperParser()
    URL_QUERY = 'http://api.springer.com/metadata/json?' \
                'q=issn:{issn} sort:date&' \
                's={s}&' \
                'p={p}&' \
                'api_key={key}'
    ITEMS_PER_PAGE = 20

    def journal_is_valid(self, journal):
        if super(ConsumerSpringer, self).journal_is_valid(journal):
            if journal.id_issn or journal.id_eissn:
                return True

        return False

    def consume_journal(self, journal):

        entries = []

        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status_to('consuming')

        start_date = self.get_start_date(cj)

        issn = journal.id_issn or journal.id_eissn
        s = 1
        if issn:
            try:
                while True:
                    query = self.URL_QUERY.format(issn=issn,
                                                  key=settings.CONSUMER_SPRINGER_KEY,
                                                  s=s,
                                                  p=self.ITEMS_PER_PAGE)

                    res = self.session.get(query).json()

                    min_date = min(list(map(self.get_publication_date,
                                            res.get('records'))))
                    if min_date < start_date:
                        entries += res.get('records')
                        break
                    else:
                        entries += res.get('records')
                        s += self.ITEMS_PER_PAGE

            except Exception:
                cj.status_to('error')
                raise

        cj.status_to('idle')

        return entries

    def get_publication_date(self, entry):
        return parse(entry.get('publicationDate')).date()


class ConsumerBiorxiv(SimpleCrawlerListMixin, Consumer):
    """BioRxiv Consumer"""

    TYPE = 'BIO'
    parser = BiorxivPaperParser()
    DOMAIN = 'http://biorxiv.org'
    LIST_PATTERN = '{domain}/archive?field_highwire_a_epubdate_value[value]&page={page}'
    DETAIL_PATTERN = {
        'name': 'a',
        'class': 'highwire-cite-linked-title',
    }

    def journal_is_valid(self, journal):
        if super(ConsumerBiorxiv, self).journal_is_valid(journal):
            if journal.id_oth == 'biorxiv':
                return True
        return False

    def consume_journal(self, journal):
        """Consume BioRxiv 'journal'"""
        # Update consumer_journal status
        cj = self.consumerjournal_set.get(journal=journal)
        cj.status_to('consuming')
        try:
            start_date = self.get_start_date(cj)
            entries = self.crawl_to_date(start_date)
        except Exception:
            cj.status_to('error')
            raise

        cj.status_to('idle')

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
                if self.base_coundown_period < settings.CONSUMER_MAX_DELAY:
                    self.base_coundown_period += 1
                else:
                    self.base_coundown_period = settings.CONSUMER_MAX_DELAY
            elif n_fet > 0:
                if self.base_coundown_period > settings.CONSUMER_MIN_DELAY:
                    self.base_coundown_period -= 1
                else:
                    self.base_coundown_period = settings.CONSUMER_MIN_DELAY
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


class ConsumerPubPeer(TimeStampedModel):

    last_consume_at = models.DateTimeField(null=True, blank=True)

    parser = PubPeerThreadParser()

    URL_QUERY = 'http://api.pubpeer.com/v1/publications/dump/'

    # API key
    API_KEY = settings.CONSUMER_PUBPEER_API_KEY

    def __init__(self, *args, **kwargs):
        super(ConsumerPubPeer, self).__init__(*args, **kwargs)
        # attach session
        self.session = requests.Session()
        adapter = adapters.HTTPAdapter(max_retries=3)
        self.session.mount('http://', adapter)

    def consume(self):
        """Retrieve PubPeer comments from API"""
        entries = []
        page = 1
        if self.last_consume_at:
            from_date = self.last_consume_at.timestamp() - 3600 * 24
        else:
            from_date = time.time() - \
                        3600 * 24 * settings.CONSUMER_PUBPEER_INIT_PAST
        while True:
            time.sleep(2)
            query = '{url}{page}?devkey={key}'.format(
                url=self.URL_QUERY,
                page=page,
                key=self.API_KEY
            )
            resp = self.session.get(query)
            entries += json.loads(resp.text)['publications']

            ds = min([float(c['date']) for d in entries for c in d['comments']])
            if ds < from_date:
                break
            page += 1

        return entries

    def populate(self):
        """Populate DB with PubPeer     comments"""
        # consume
        entries = self.consume()

        # save to database
        count = 0
        for entry in entries:
            item = self.parser.parse(entry)
            if item['pubpeer']['doi']:
                try:
                    thread = self.add_or_update_entry.delay(item)
                except RuntimeError:
                    thread = None
                    pass
                if thread:
                    count += 1
        self.last_consume_at = timezone.now()
        self.save()
        return count

    @app.task(filter=task_method)
    def add_or_update_entry(self, entry):

        ppm = PubPeerManager()
        thread = ppm.add_or_update_entry(entry)

        return thread
