from __future__ import print_function
import csv
import logging
from progressbar import ProgressBar, Percentage, Bar, ETA
from django.db.models import Q
from paperstream.library.models import Publisher, Journal
from paperstream.library.forms import JournalFormFillBlanks, PublisherForm
from paperstream.consumers.models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier
from paperstream.consumers.constants import CONSUMER_TYPE

logger = logging.getLogger(__name__)

def populate_journal(csv_file, print_to=None):

    # counting number of row
    with open(csv_file, 'rU') as rows:
        reader = csv.DictReader(rows, delimiter=';')
        trows = sum(1 for row in reader)

    # Init
    records_added = 0
    errors = []
    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=trows, redirect_stderr=True).start()

    with open(csv_file, 'rU') as rows:
        # Generate a dict per row, with the first CSV row being the keys.
        for i, row in enumerate(csv.DictReader(rows, delimiter=";")):
            try:
                try:
                    journal = Journal.objects.get(
                        Q(id_issn=row['id_issn']) |
                        Q(id_eissn=row['id_eissn']) |
                        Q(id_oth=row['id_oth']) |
                        Q(id_arx=row['id_arx']))
                except Journal.DoesNotExist:
                    journal = None

                # JournalFormFillUp is used here because we got multiple source of
                # journals that have complementary information
                form = JournalFormFillBlanks(row, instance=journal)
                if form.is_valid():
                    journal = form.save()
                    # Save as trusted
                    journal.is_trusted = True
                    journal.save()
                    records_added += 1
                else:
                    errors.append(form.errors)
                    logger.warning('journal failed to register: '
                                   '{0}'.format(form.errors))
            except KeyError as e:
                logger.warning('journal row error: {0}'.format(row))
                pass

                # update progress bar
            pbar.update(i)
        # close progress bar
        pbar.finish()

    return records_added, errors


def populate_publisher(csv_file, print_to=None):

    records_added = 0
    errors = []

    # counting the number of lines. Useful for progressbar
    # counting number of row
    with open(csv_file, 'r') as rows:
        reader = csv.DictReader(rows, delimiter=';')
        trows = sum(1 for row in reader)

    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=trows, redirect_stderr=True).start()

    with open(csv_file, 'r') as rows:
        for i, row in enumerate(csv.DictReader(rows, delimiter=";")):
            try:
                publisher = Publisher.objects.get(name=row['name'])
            except Publisher.DoesNotExist:
                publisher = None

            form = PublisherForm(row, instance=publisher)
            if form.has_changed():
                if form.is_valid():
                    form.save()
                    records_added += 1
                else:
                    errors.append(form.errors)

            # update progress bar
            pbar.update(i)
        # close progress bar
        pbar.finish()

    return records_added, errors


def populate_consumer(type_, name, csv_file=None, print_to=None):

    if type_ == 'PUB':
        consumer = ConsumerPubmed.objects.create(name=name)
        records_added, errors = populate_consumer_from_file(consumer,
                                                            type_,
                                                            csv_file,
                                                            print_to=print_to)
    elif type_ == 'ARX':
        consumer = ConsumerArxiv.objects.create(name=name)
        records_added, errors = populate_consumer_from_db(consumer,
                                                          type_,
                                                          print_to=print_to)
    elif type_ == 'ELS':
        consumer = ConsumerElsevier.objects.create(name=name)
        records_added, errors = populate_consumer_from_db(consumer,
                                                          type_,
                                                          print_to=print_to)
    else:
        raise ValueError('unknown type {0}'.format(type_))

    # activate all entry in journals
    consumer.activate_all()

    return records_added, errors

def populate_consumer_from_db(consumer, type_, print_to=None):

    # Init
    records_added = 0
    errors = []

    # Get journal from db
    publisher_name = dict(CONSUMER_TYPE)[type_]
    journals = Journal.objects.filter(publisher__name=publisher_name)

    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=journals.count(), redirect_stderr=True).start()

    for i, journal in enumerate(journals):
        try:
            consumer.add_journal(journal)
            records_added += 1
        except Exception as e:
            logging.warning('journal failed to link to consumer {0}'
                            ''.format(journal.title))
            pass
        # update progress bar
        pbar.update(i)
    # close progress bar
    pbar.finish()

    return records_added, errors

def populate_consumer_from_file(consumer, type_, csv_file, print_to=None):

    # Init
    records_added = 0
    errors = []

    # counting number of row
    with open(csv_file, 'r') as rows:
        reader = csv.DictReader(rows, delimiter=';')
        trows = sum(1 for row in reader)

    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=trows, redirect_stderr=True).start()

    with open(csv_file, 'r') as rows:
        # Generate a dict per row, with the first CSV row being the keys.
        for i, row in enumerate(csv.DictReader(rows, delimiter=";")):
            try:
                if type_ == 'PUB':
                    journal = Journal.objects.get(id_issn=row['id_issn'])
                elif type_ == 'ARX':
                    journal = Journal.objects.get(id_arx=row['id_arx'])
                else:
                    journal = None

                consumer.add_journal(journal)

                records_added += 1

            except Exception as e:
                logging.warning('journal failed to link to consumer row {0}'.format(row))
                pass
            # update progress bar
            pbar.update(i)

        # close progress bar
        pbar.finish()

    return records_added, errors