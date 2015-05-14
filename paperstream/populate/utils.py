import csv
import sys
import tempfile
from progressbar import ProgressBar, Percentage, Bar, ETA
from library.models import Publisher, Journal
from library.forms import JournalFormFillUp, PublisherForm
from django.db.models import Q


def populate_journal(csv_file, print_to=None):

    # counting number of row
    with open(csv_file, 'r') as rows:
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

    with open(csv_file, 'r') as rows:
        # Generate a dict per row, with the first CSV row being the keys.
        for i, row in enumerate(csv.DictReader(rows, delimiter=";")):
            try:
                journal = Journal.objects.get(Q(id_issn=row['id_issn']) |
                                    Q(id_eissn=row['id_eissn']) |
                                    Q(id_oth=row['id_oth']) |
                                    Q(id_arx=row['id_arx']))
            except Journal.DoesNotExist:
                journal = None

            # JournalFormFillUp is used here because we got multiple source of
            # journals that have complementary information
            form = JournalFormFillUp(row, instance=journal)
            if form.is_valid():
                journal = form.save()
                # Save as trusted
                journal.is_trusted = True
                journal.save()
                records_added += 1
            else:
                errors.append(form.errors)

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
