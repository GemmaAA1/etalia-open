import csv
import tempfile
from django import forms
from paperstream.library.models import Journal, Publisher
from progressbar import ProgressBar, Percentage, Bar, ETA


def populate_publisher(rows):

    records_added = 0
    errors = []

    for row in csv.DictReader(rows, delimiter=","):

        form = PublisherForm(row)
        if form.is_valid():
            form.save()
            records_added += 1
        else:
            errors.append(form.errors)

    return records_added, errors


def populate_journal(rows):

    records_added = 0
    errors = []
    # Generate a dict per row, with the first CSV row being the keys.
    for row in csv.DictReader(rows, delimiter=";"):
        # Bind the row data to the JournalThomsonReutersForm.
        form = JournalForm(row)
        # Check to see if the row data is valid.
        if form.is_valid():
            form.save()
            records_added += 1
        else:
            errors.append(form.errors)
    return records_added, errors


def preprocess_thomsonreuters_file(in_csv_file):

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
    # header of output file
    fieldnames = ['id_issn', 'id_eissn', 'title', 'period', 'publisher']
    # output file
    out_csv_file = tempfile.NamedTemporaryFile().name+'.csv'
    # counting the number of lines. Useful for progressbar
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        trows = sum(1 for row in reader)


    with open(out_csv_file, 'w') as csv_processed:
        writer = csv.DictWriter(csv_processed, fieldnames=fieldnames)

        with open(input_file, 'r') as in_csv_file:
            # loop through input file
            for i, row in enumerate(csv.DictReader(in_csv_file, delimiter=";")):

                # init progress bar
                pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                                   maxval=trows).start()

                # clean period
                period = row['period']
                if period.lower() in PUBLISH_PERIOD_MATCH.keys():
                    period = PUBLISH_PERIOD_MATCH[period.lower()]
                else:
                    period = ''

                # clean publisher
                # publisher is csv files is long string with the name and the address
                # that needs to be fetched
                publisher = row['publisher']
                publishers = Publisher.objects.all()
                publisher_match = [pub for pub in publishers if pub.name
                                   in publisher]
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

                # Write preprocess row to out file
                writer.writerow({'id_issn': row['id_issn'],
                                 'id_eissn': row['id_eissn'],
                                 'title': title,
                                 'period': period,
                                 'publisher': publisher})
                # update progress bar
                pbar.update(i)

    return out_csv_file


if __name__ == '__main__':

    # Populate publisher
    input_file = 'data/publishers/publisher_list.csv'
    with open(input_file) as csv_file:
        records_added, errors = populate_publisher(csv_file)

    # Populate journals
    input_file = 'data/journals/20150510_thomsonreuters.csv'
    # first, preprocess to match pipe
    out_file = preprocess_thomsonreuters_file(input_file)
    # then, populate
    with open(out_file) as csv_file:
        records_added, errors = populate_journal(csv_file)