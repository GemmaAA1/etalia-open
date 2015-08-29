# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.test import TestCase
from django.utils import timezone

from paperstream.library.models import Journal, Paper


class LibDataTestCase(TestCase):

    def setUp(self):
        self.journal = Journal.objects.create(title='Journal test')
        self.paper = Paper.objects.create(
            title='Bla bla.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date())
        self.paper2 = Paper.objects.create(
            title='Blo blo.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date())
        self.paper3 = Paper.objects.create(
            title='Bli bli.',
            abstract='Hi. Hi, <p>hi</p> {mu}\n',
            journal=self.journal,
            date_ep=timezone.now().date())