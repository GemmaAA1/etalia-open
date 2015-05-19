import os
from django.conf import settings
from populate.utils import populate_publisher, populate_journal, \
    populate_consumer
from django.core.management.base import BaseCommand
from consumers.constants import CONSUMER_TYPE
from ...constants import PUBLISHER_OPTIONS, JOURNAL_OPTIONS, CONSUMER_OPTIONS

class Command(BaseCommand):
    """

    """

    help = 'Pre-Populate paperstream Journal (populate jounral), Publisher ' \
           '(populate publisher), Consumer (populate consumer) or all ' \
           '(populate all)'

    def add_arguments(self, parser):
        parser.add_argument('what', nargs=1, type=str)
        parser.add_argument('source', nargs='?', type=str, default='all')
        parser.add_argument('-n', '--name', action='store', dest='name',
                            help='Consumer name', type=str)

    def handle(self, *args, **options):
        pop_what = options['what'][0]
        source = options['source']
        name = options['name']

        self.stdout.write(
            'Populating {what}'.format(what=pop_what.capitalize()))

        # populate publisher
        if pop_what == 'publisher':
            if source == 'all':
                # Populate publisher
                for sourcefile in PUBLISHER_OPTIONS:
                    self.stdout.write(
                        '... from {0}'.format(sourcefile['source'].capitalize()))
                    records_added, errors = populate_publisher(
                        sourcefile['file_path'],
                        print_to=self.stderr)
                    # print(errors)
            else:
                msg = 'argument {opt} unknown, possible choice ' \
                      'is: {ch}'.format(
                    opt=options['source'],
                    ch=' '.join([s['source'] for s in PUBLISHER_OPTIONS]))
                raise ValueError(msg)

        # populate journals
        elif pop_what == 'journal':

            source_options = [s['source'] for s in JOURNAL_OPTIONS]
            print_source_options = source_options + ['all']

            # Populate journals
            if source == 'all':
                for sourcefile in JOURNAL_OPTIONS:
                    self.stdout.write(
                        '... from {0}'.format(
                            sourcefile['source'].capitalize()))
                    records_added, errors = populate_journal(
                        sourcefile['file_path'],
                        print_to=self.stderr)
                    # print(errors)
            elif source in source_options:
                sourcefile = [so for so in JOURNAL_OPTIONS if
                              so['source'] == source][0]
                self.stdout.write(
                    '... from {0}'.format(
                        sourcefile['source'].capitalize()))
                records_added, errors = populate_journal(
                        sourcefile['file_path'],
                        print_to=self.stderr)
            else:
                msg = 'argument {opt} unknown, possible choice ' \
                      'is: {ch}'.format(
                    opt=source,
                    ch=' '.join(print_source_options))
                raise ValueError(msg)

        elif pop_what == 'consumer':

            source_options = [h_read.lower() for code, h_read in CONSUMER_TYPE]
            source_dict = dict([(cons[1].lower(), cons[0])
                                for cons in CONSUMER_TYPE])
            print_source_options = source_options

            if source in source_options:
                self.stdout.write(
                    '... from {0}'.format(source.capitalize()))
                try:
                    sourcefile = [so for so in CONSUMER_OPTIONS if
                          so['source'] == source][0]
                except IndexError:
                    sourcefile = None
                # if file exist, then populate with file
                if sourcefile:
                    file = sourcefile['file_path']
                else:
                    file = None

                populate_consumer(source_dict[source], name, file,
                                  print_to=self.stderr)
            else:
                msg = 'argument {opt} unknown, possible choice is: {ch}'.format(
                    opt=source,
                    ch=' '.join(print_source_options))
                raise ValueError(msg)

        else:
            raise ValueError('Unknown argument: {0}'.format(pop_what))




