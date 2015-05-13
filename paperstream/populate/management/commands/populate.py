import os
from django.conf import settings
from populate.utils import populate_publisher, populate_journal
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Pre-Populate paperstream Journal (populate jounral), Publisher ' \
           '(populate publisher), Consumer (populate consumer) or all ' \
           '(populate all)'

    publisher_list_file = [
        {'source': 'Home',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/publishers/publisher_list.csv')},
    ]

    journal_list_file = [
        {'source': 'Thomson Reuters Master List',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_thomsonreuters_cleaned.csv')},
        {'source': 'MedLine List',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_medline_cleaned.csv')},
        {'source': 'Arxiv home build list',
         'file_path': os.path.join(str(settings.STATIC_ROOT),
                     'populate/journals/20150510_arxiv_cleaned.csv')},
    ]

    def add_arguments(self, parser):
        parser.add_argument('populate_what', nargs='+', type=str)

    def handle(self, *args, **options):
        for pop_what in options['populate_what']:

            # populate publisher
            if pop_what == 'publisher':
                self.stdout.write(
                    'Populating {what}'.format(what=pop_what.capitalize()))
                # Populate publisher
                for sourcefile in self.publisher_list_file:
                    self.stdout.write(
                        '... from {0}'.format(sourcefile['source']))
                    records_added, errors = populate_publisher(
                        sourcefile['file_path'],
                        print_to=self.stderr)
                    # print(errors)

            # populate journals
            if pop_what == 'journal':
                self.stdout.write(
                    'Populating {what}'.format(what=pop_what.capitalize()))
                # Populate journals
                for sourcefile in self.journal_list_file:
                    self.stdout.write(
                        '... from {0}'.format(sourcefile['source']))
                    records_added, errors = populate_journal(
                        sourcefile['file_path'],
                        print_to=self.stderr)
                    # print(errors)




