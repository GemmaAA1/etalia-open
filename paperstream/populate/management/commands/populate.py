import os
from django.conf import settings
from populate.utils import preprocess_thomsonreuters_file, populate_publisher, \
    populate_journal
from django.core.management.base import BaseCommand, CommandError



class Command(BaseCommand):
    help = 'Pre-Populate paperstream Journal (populate jounral), Publisher ' \
           '(populate publisher), Consumer (populate consumer) or all ' \
           '(populate all)'

    def add_arguments(self, parser):
        parser.add_argument('populate_what', nargs='+', type=str)

    def handle(self, *args, **options):
        for pop_what in options['populate_what']:

            if pop_what == 'publisher':
                self.stdout.write('Populating {what}'.format(what='Publishers'))
                # Populate publisher
                input_file = settings.STATIC_ROOT + \
                             '/populate/publishers/publisher_list.csv'
                records_added, errors = populate_publisher(
                    input_file,
                    print_to=self.stderr)

            if pop_what == 'journal':
                self.stdout.write('Populating {what}'.format(what='Journals'))
                # Populate journals

                # Master journal list from Thomson and Reuters (WebOfKnowledge)
                self.stdout.write('...from {0}'.format('Thomson Reuters list'))
                self.stdout.write('...preprocessing')
                input_file = settings.STATIC_ROOT + \
                             '/populate/journals/20150510_thomsonreuters.csv'
                # input_file = settings.STATIC_ROOT + \
                #              '/populate/journals/test.csv'
                # first, preprocess to match pipe
                output_file = preprocess_thomsonreuters_file(
                    input_file,
                    print_to=self.stderr)
                # then, populate
                self.stdout.write('...populating')
                records_added, errors = populate_journal(output_file)
                os.remove(output_file)




