#!/usr/bin/env python
import csv
import paperstream.apps.library.utils as libu
import paperstream.apps.library.models as libm
import paperstream.apps.consumers.models as conm
from progressbar import ProgressBar, Percentage, Bar, ETA


def pop_publishers():
    print("Starting db population of publishers...")
    libu.add_publisher(name='ELSEVIER',
                       url='http://www.elsevier.com')
    libu.add_publisher(name='ARXIV',
                       url='http://arxiv.org')
    libu.add_publisher(name='WILEY',
                       url='http://onlinelibrary.wiley.com/')
    libu.add_publisher(name='PUBLIC LIBRARY SCIENCE',
                       url='http://www.plos.org/')

    # Add other known publisher here...

def pop_consumers():

    #Elsevier
    pub = libm.Publisher.objects.get(name='ELSEVIER')
    cons = conm.Consumer.objects.get_or_create(name='ELSEVIER')[0]
    cons.type = 'API'
    cons.publisher = pub
    cons.conf_file = './static/config/consumers/elsevier.json',
    cons.is_active = True
    cons.save()

    # Arxiv
    cons = conm.Consumer.objects.get_or_create(name='ARXIV')[0]
    cons.type = 'API'
    cons.publisher = None
    cons.conf_file = './static/config/consumers/arxiv.json'
    cons.is_active = True
    cons.save()

    # Pubmed
    cons = conm.Consumer.objects.get_or_create(name='PUBMED')[0]
    cons.type = 'API'
    cons.publisher = None
    cons.conf_file = './static/config/consumers/pubmed.json'
    cons.is_active = True
    cons.save()


def pop_journals():
    print("Starting db population of journals...")

    # Publisher
    pubs = libm.Publisher.objects.all()

    # Loop on master journal list data file
    print("from Thomson & Reuters...")
    input_file = './static/fixtures/150326_mjl_consolidated.csv'
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        trows = sum(1 for row in reader)

    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=trows).start()
        cou = 0
        for row in reader:
            if row['issn']:
                j, j_new = libu.add_known_journal(issn=row['issn'])

                if j_new:
                    j.e_issn = row.get('e_issn', '')
                    j.title = row.get('title', '')
                    j.short_title = row.get('short_title', '')
                    j.period = row.get('period', '')

                    p_candidate = [p for p in pubs if p.name in row['publisher']]
                    if len(p_candidate) == 1:
                        j.publisher = p_candidate[0]
                    else:
                        j.publisher = None

                    j.save()

                # consume journal
                if j.publisher and \
                   conm.Consumer.objects.filter(publisher=j.publisher).exists():

                    c = conm.Consumer.objects.get(publisher=j.publisher)
                    conm.ConsumerJournal.objects.get_or_create(
                        journal=j,
                        consumer=c,
                        is_active=False)
                # update progress bar
            cou += 1
            pbar.update(cou)
        pbar.finish()


    # Arxiv - Arxiv is split in different journal for each category
    print("from Arxiv...")
    input_file = './static/fixtures/arxiv_ext_id_2015.csv'
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        trows = sum(1 for row in reader)
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=trows).start()
        cou = 0
        for row in reader:
            j, j_new = libu.add_known_journal(ext_id=row['ext_id'])

            if j_new:
                j.title = row.get('title', '')
                j.url = row.get('url', '')
                j.publisher = libm.Publisher.objects.get(name='ARXIV')
                j.save()

            # consume journal
            if conm.Consumer.objects.filter(name='ARXIV').exists():
                c = conm.Consumer.objects.get(name='ARXIV')
                conm.ConsumerJournal.objects.get_or_create(
                    journal=j,
                    consumer=c,
                    is_active=False)

            cou += 1
            pbar.update(cou)
        pbar.finish()


    # Pubmed -
    print("from Pubmed...")
    input_file = './static/fixtures/150328_jl_medline.csv'
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        trows = sum(1 for row in reader)
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=';')
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=trows).start()
        cou = 0
        for row in reader:
            j, j_new = libu.add_known_journal(issn=row.get('issn', ''))
            if j_new:
                j.e_issn = row.get('e_issn', '')
                j.title = row.get('title', '')
                j.short_title = row.get('short_title', '')
                j.publisher = None
                j.save()

            # consume journal
            if conm.Consumer.objects.filter(name='PUBMED').exists():
                c = conm.Consumer.objects.get(name='PUBMED')
                conm.ConsumerJournal.objects.get_or_create(
                    journal=j,
                    consumer=c,
                    is_active=False)

            cou += 1
            pbar.update(cou)
        pbar.finish()


if __name__ == '__main__':

    # Order matters...
    pop_publishers()
    pop_consumers()
    pop_journals()