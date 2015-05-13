import csv
import sys
import tempfile
from progressbar import ProgressBar, Percentage, Bar, ETA
from library.models import Publisher, Journal
from library.forms import JournalForm, PublisherForm
from django.db.models import Q

def preprocess_thomsonreuters_file(in_csv_file, print_to=None):

    # Init
    # To deal with disparity in dataset and match model period choices
    PUBLISH_PERIOD_MATCH = {'annual': 'ANN',
                            'semiannual': 'SEM',
                            'semi-annual': 'SEM',
                            'quarterly': 'QUA',
                            'tri-annual': 'TRI',
                            'triannual': 'TRI',
                            'bi-monthly': 'BIM',
                            'bimonthly': 'BIM',
                            'monthly': 'MON',
                            'irregular': 'IRR'}

    # Retrieve all publishers
    publishers = Publisher.objects.all()

    # header of output file
    fieldnames = ['id_issn', 'id_eissn', 'title', 'short_title',
                  'period', 'publisher']
    # output file
    out_csv_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.csv')
    # counting the number of lines. Useful for progressbar
    with open(in_csv_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        trows = sum(1 for row in reader)

    # Process input file
    # Open output file
    # with open(out_csv_file, 'w') as csv_processed:
    writer = csv.DictWriter(out_csv_file, fieldnames=fieldnames,
                            delimiter=';')
    writer.writeheader()

    # open input file
    with open(in_csv_file, 'r') as in_csv_file:
        # init progress bar
        # for redirecting progress bar to proper stderr
        if print_to:
            print('', file=print_to)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=trows, redirect_stderr=True).start()

        # loop through input file
        for i, row in enumerate(csv.DictReader(in_csv_file, delimiter=";")):

            # clean period
            period = row['period']
            if period.lower() in PUBLISH_PERIOD_MATCH.keys():
                period = PUBLISH_PERIOD_MATCH[period.lower()]
            else:
                period = ''

            # clean publisher
            # publisher is csv files is long string with the name and the address
            # that needs to be fetched
            publisher = row['publisher'].lower()

            publisher_match = [pub for pub in publishers if pub.name.lower() in
                               publisher]
            # if multiple match, reject
            if len(publisher_match) == 1:
                publisher = publisher_match[0].name
            else:
                publisher = None

            # clean title
            title = row['title']
            # if title all upper or lower, capitalized
            if title.isupper() or title.islower():
                title = title.title()

            # clean short_title
            short_title = row['short_title']
            # if title all upper or lower, capitalized
            if short_title.isupper() or short_title.islower():
                short_title = short_title.title()

            # Write preprocess row to out file
            writer.writerow({'id_issn': row['id_issn'],
                             'id_eissn': row['id_eissn'],
                             'title': title,
                             'short_title': short_title,
                             'period': period,
                             'publisher': publisher})

            # update progress bar
            pbar.update(i)
        # close progress bar
        pbar.finish()
    output_filename = out_csv_file.name
    out_csv_file.close()
    return output_filename


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

            # Retrieve journal data if ids are found
            try:
                journal = Journal.objects.get(Q(id_issn=row['id_issn']) |
                                    Q(id_eissn=row['id_eissn']) |
                                    Q(id_oth=row['id_oth']) |
                                    Q(id_arx=row['id_arx']))
            except Journal.DoesNotExist:
                journal = None

            # Complete row field if empty
            for field in JournalForm.Meta.fields:
                if not row[field] and journal:
                    row[field] = getattr(journal, field)

            form = JournalForm(row, instance=journal)
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
