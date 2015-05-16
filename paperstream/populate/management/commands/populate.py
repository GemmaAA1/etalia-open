import os
from django.conf import settings
from populate.utils import populate_publisher, populate_journal
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    """

    """

    help = 'Pre-Populate paperstream Journal (populate jounral), Publisher ' \
           '(populate publisher), Consumer (populate consumer) or all ' \
           '(populate all)'

    publisher_options = [
        {'source': 'all',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/publishers/publisher_list.csv')},
    ]

    journal_options = [
        {'source': 'thomson',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_thomsonreuters_cleaned.csv')},
        {'source': 'pubmed',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_medline_cleaned.csv')},
        {'source': 'arxiv',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_arxiv_cleaned.csv')},
    ]

    def add_arguments(self, parser):
        parser.add_argument('what', nargs=1, type=str, default='all')
        parser.add_argument('source', nargs=1, type=str, default='all')

    def handle(self, *args, **options):
        pop_what = options['what'][0]
        source = options['source'][0]

        # populate publisher
        if pop_what == 'publisher':
            if source == 'all':
                self.stdout.write(
                    'Populating {what}'.format(what=pop_what.capitalize()))
                # Populate publisher
                for sourcefile in self.publisher_options:
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
                    ch=' '.join([s['source'] for s in self.publisher_options]))
                raise ValueError(msg)

        # populate journals
        if pop_what == 'journal':
            self.stdout.write(
                'Populating {what}'.format(what=pop_what.capitalize()))

            source_options = [s['source'] for s in self.journal_options]
            print_source_options = source_options + ['all']

            # Populate journals
            if source == 'all':
                for sourcefile in self.journal_options:
                    self.stdout.write(
                        '... from {0}'.format(
                            sourcefile['source'].capitalize()))
                    records_added, errors = populate_journal(
                        sourcefile['file_path'],
                        print_to=self.stderr)
                    # print(errors)
            elif source in source_options:
                sourcefile = [so for so in self.journal_options if
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




