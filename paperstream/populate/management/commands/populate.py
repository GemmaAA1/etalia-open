# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from paperstream.populate.utils import populate_publisher, populate_journal, \
    populate_consumer
from django.core.management.base import BaseCommand
from paperstream.consumers.constants import CONSUMER_TYPE
from ...constants import PUBLISHER_OPTIONS, JOURNAL_OPTIONS, CONSUMER_OPTIONS

class Command(BaseCommand):
    """
    Pre-Populate paperstream models: Three models <what> can be populated:
       - Journal (populate journal)\n
       - Publisher (populate publisher) \n
       - Consumer (populate consumer)\n
    The <source> for pre-population must be specified and defined in constant.py.
    For Journal and Publisher, <source>=all will populate models with all
    possible sources as defined in constant.py
    Note: journal associated to consumers are all active when populated

    Ex:
    >> ./manage populate journal thomson
    will create Journal instances based on the thomson list file as defined in constant.py

    >> ./manage populate consumer pubmed --name pubmed1
    will create one PubmedConsumer based on pubmed list as defined in constant.py
    """

    help = 'Pre-Populate paperstream models'

    def add_arguments(self, parser):
        parser.add_argument('what', nargs=1, type=str,
                            help='specify the model (journal, publisher, consumer or all)')
        parser.add_argument('source', nargs='?', type=str, default='all',
                            help='specify source name (all or as defined in constant.py')
        parser.add_argument('-n', '--name', action='store', dest='name',
                            help='specify consumer name (mandatory if populating consumer)', type=str)
        parser.add_argument('-l', '--local', action='store', dest='local',
                            help='local flag', type=bool, default=False)

    def handle(self, *args, **options):
        pop_what = options['what'][0]
        source = options['source']
        name = options['name']
        local = options['local']

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
                          so['source'] == source and so['local'] == local][0]
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




