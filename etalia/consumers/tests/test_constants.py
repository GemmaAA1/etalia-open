# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.test import TestCase
from etalia.library.constants import PAPER_TYPE
from ..constants import PUBMED_PT

class ConstantTest(TestCase):

    def test_pubmed_type_subset_paper_type(self):
        pubmed_type = set(dict(PUBMED_PT).values())
        paper_type = set(dict(PAPER_TYPE).keys())
        self.assertTrue(paper_type.issuperset(pubmed_type))

