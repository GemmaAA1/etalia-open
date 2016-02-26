# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import, print_function

import csv
import json
import logging
from progressbar import ProgressBar, Percentage, Bar, ETA
from django.db.models import Q
from paperstream.library.models import Publisher, Journal
from paperstream.library.forms import JournalFormFillBlanks, PublisherForm
from paperstream.consumers.models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier
from paperstream.consumers.constants import CONSUMER_TYPE

logger = logging.getLogger(__name__)


def populate_journal(json_file, print_to=None):

    # counting number of row
    with open(json_file, 'r') as f:
        data = json.load(f)
    nb_elements = len(data)

    # Init
    records_added = 0
    errors = []
    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=nb_elements, redirect_stderr=True).start()

    for i, el in enumerate(data):
        try:
            try:
                journal = Journal.objects.get(
                    Q(id_issn=el['id_issn']) |
                    Q(id_eissn=el['id_eissn']) |
                    Q(id_oth=el['id_oth']) |
                    Q(id_arx=el['id_arx']))
            except Journal.DoesNotExist:
                journal = None

            # JournalFormFillUp is used here because we got multiple source of
            # journals that have complementary information
            form = JournalFormFillBlanks(el, instance=journal)
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
        except Exception as e:
            logger.warning('journal row error: {0}'.format(el))
            pass

            # update progress bar
        pbar.update(i)
    # close progress bar
    pbar.finish()

    return records_added, errors


def populate_publisher(json_file, print_to=None):

    records_added = 0
    errors = []

    with open(json_file, 'r') as f:
        data = json.load(f)
    nb_elements = len(data)

    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=nb_elements, redirect_stderr=True).start()

    for i, el in enumerate(data):
        try:
            publisher = Publisher.objects.get(name=el['name'])
        except Publisher.DoesNotExist:
            publisher = None

        form = PublisherForm(el, instance=publisher)
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


def populate_consumer(type_, name, json_file=None, print_to=None):

    if type_ == 'PUB':
        consumer, _ = ConsumerPubmed.objects.get_or_create(name=name)
        records_added, errors = populate_consumer_from_file(consumer,
                                                            type_,
                                                            json_file,
                                                            print_to=print_to)
    elif type_ == 'ARX':
        consumer, _ = ConsumerArxiv.objects.get_or_create(name=name)
        records_added, errors = populate_consumer_from_db(consumer,
                                                          type_,
                                                          print_to=print_to)
    elif type_ == 'ELS':
        consumer, _ = ConsumerElsevier.objects.get_or_create(name=name)
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


def populate_consumer_from_file(consumer, type_, json_file, print_to=None):

    # Init
    records_added = 0
    errors = []

    with open(json_file, 'r') as f:
        data = json.load(f)
    nb_elements = len(data)

    # init progress bar
    # for redirecting progress bar to proper stderr
    if print_to:
        print('', file=print_to)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                       maxval=nb_elements, redirect_stderr=True).start()

    for i, el in enumerate(data):
        try:
            if type_ == 'PUB':
                try:
                    journal = Journal.objects.get(id_issn=el['id_issn'])
                except Journal.DoesNotExist:
                    journal = Journal.objects.get(id_eissn=el['id_eissn'])
            elif type_ == 'ARX':
                journal = Journal.objects.get(id_arx=el['id_arx'])
            else:
                journal = None

            consumer.add_journal(journal)

            records_added += 1

        except Exception as e:
            logging.warning('journal failed to link to consumer row {0}'.format(el))
            pass
        # update progress bar
        pbar.update(i)

    # close progress bar
    pbar.finish()

    return records_added, errors