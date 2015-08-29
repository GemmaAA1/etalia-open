# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.test import TestCase
from dateutil.parser import parse

from ..parsers import ParserPubmed, ParserArxiv, ParserElsevier

# Not sure what the best way to unit test parser.
# If test is ran on multiple entries, testing code might really
# look like parser code if we do not manually enter all the entry fields


class ParserPubmedTest(TestCase):

    def setUp(self):
        with open('paperstream/consumers/tests/pubmed_sample.json') as file:
            self.entries = json.load(file)

    def test_parse_journal(self):
        # test only first record
        entry = self.entries[0]
        parser = ParserPubmed()
        journal = parser.parse_journal(entry)
        self.assertEqual(journal['title'], 'Neurosurgery')
        self.assertEqual(journal['id_issn'], '0148-396X')
        self.assertEqual(journal['id_eissn'], '1524-4040')

    def test_parse_paper(self):
        # test only first record
        entry = self.entries[0]
        parser = ParserPubmed()
        paper = parser.parse_paper(entry)
        self.assertEqual(paper['id_doi'], '10.1227/NEU.0000000000000714')
        self.assertEqual(paper['id_pii'], '00006123-201506000-00024')
        self.assertEqual(paper['abstract'][:len('BACKGROUND:')], 'BACKGROUND:')
        self.assertEqual(paper['type'], 'JOU')
        self.assertEqual(paper['id_pmi'], '25988929')
        self.assertEqual(paper['publish_status'], 'ppublish')
        self.assertEqual(paper['volume'], '76')
        self.assertEqual(paper['issue'], '6')
        self.assertEqual(paper['page'], '756-765')
        self.assertEqual(paper['date_ep'], None)
        self.assertEqual(paper['date_pp'], parse('2015 Jun'))
        self.assertEqual(paper['date_lr'], parse('20150520'))
        self.assertEqual(paper['url'], 'http://doi.org/10.1227/NEU.0000000000000714')
        self.assertEqual(paper['language'], 'ENG')
        self.assertEqual(paper['source'], '')

    def test_parse_authors(self):
        entry = self.entries[0]
        parser = ParserPubmed()
        authors = parser.parse_authors(entry)
        self.assertEqual(authors[0]['first_name'], 'Srivatsan')
        self.assertEqual(authors[0]['last_name'], 'Pallavaram')
        self.assertEqual(authors[1]['first_name'], 'Pierre-Francois')
        self.assertEqual(authors[1]['last_name'], "D'Haese")
        self.assertEqual(authors[-1]['first_name'], 'Joseph S')
        self.assertEqual(authors[-1]['last_name'], 'Neimat')

    def test_parse_corp_author(self):
        pass

    def test_parse_all(self):
        records = []
        parser = ParserPubmed()
        for entry in self.entries:
            records.append(parser.parse(entry))
        self.assertEqual(len(records), len(self.entries))


class ParserElsevierTest(TestCase):

    def setUp(self):
        with open('paperstream/consumers/tests/elsevier_sample.json') as file:
            self.entries = json.load(file)

    def test_parse_journal(self):
        # test only first record
        entry = self.entries[0]
        parser = ParserElsevier()
        journal = parser.parse_journal(entry)
        self.assertEqual(journal['title'], 'European Journal of Medical Genetics')
        self.assertEqual(journal['id_issn'], '1769-7212')
        self.assertEqual(journal['id_eissn'], '')

    def test_parse_paper(self):
        # test only first record
        entry = self.entries[0]
        parser = ParserElsevier()
        paper = parser.parse_paper(entry)
        self.assertEqual(paper['id_doi'], '10.1016/j.ejmg.2015.05.002')
        self.assertEqual(paper['id_pii'], 'S1769721215000920')
        self.assertEqual(paper['abstract'][:len('AbstractWilliams')],
                         'AbstractWilliams')
        self.assertEqual(paper['type'], 'JOU')
        self.assertEqual(paper['id_pmi'], '')
        self.assertEqual(paper['publish_status'], 'aheadofprint')
        self.assertEqual(paper['volume'], '')
        self.assertEqual(paper['issue'], '0')
        self.assertEqual(paper['page'], '')
        self.assertEqual(paper['date_ep'], parse('20 May 2015'))
        self.assertEqual(paper['date_pp'], None)
        self.assertEqual(paper['date_lr'], None)
        self.assertEqual(paper['url'],
            'http://www.sciencedirect.com/science/article/pii/S1769721215000920X')
        self.assertEqual(paper['language'], 'ENG')
        self.assertEqual(paper['source'], '')

    def test_parse_authors(self):
        entry = self.entries[0]
        parser = ParserElsevier()
        authors = parser.parse_authors(entry)
        self.assertEqual(authors[0]['first_name'], 'Kimiko')
        self.assertEqual(authors[0]['last_name'], 'Ueda')
        self.assertEqual(authors[1]['first_name'], 'Junji')
        self.assertEqual(authors[1]['last_name'], "Yamada")
        self.assertEqual(authors[-1]['first_name'], 'Nobuhiko')
        self.assertEqual(authors[-1]['last_name'], 'Okamoto')

    def test_parse_corp_author(self):
        pass

    def test_parse_all(self):
        records = []
        parser = ParserElsevier()
        for entry in self.entries:
            records.append(parser.parse(entry))
        self.assertEqual(len(records), len(self.entries))


class ParserArxivTest(TestCase):
# TODO: Test ConsumerArxiv

    def test_parserarxiv_journal(self):
        pass

    def test_parserarxiv_paper(self):
        pass

    def test_parserarxiv_authors(self):
        pass

    def test_parserarxiv_corp_author(self):
        pass