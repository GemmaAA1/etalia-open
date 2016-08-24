# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
import requests
import time
import datetime
from django.db import models
from django.conf import settings
from django.contrib.auth import get_user_model
from etalia.consumers.utils import Consolidate, PaperManager

from etalia.core.parsers import PaperParser
from etalia.threads.forms import PubPeerCommentForm, PubPeerForm, ThreadForm
from etalia.core.models import TimeStampedModel
from etalia.library.models import Paper

from .parsers import PubPeerThreadParser
from .models import Thread, PubPeerComment, PubPeer

User = get_user_model()


class PubPeerConsumer(TimeStampedModel):

    last_consume_at = models.DateTimeField(null=True, blank=True)

    parser = PubPeerThreadParser()

    URL_QUERY = 'http://api.pubpeer.com/v1/publications/dump/'

    # API key
    API_KEY = settings.CONSUMER_PUBPEER_API_KEY

    def consume(self):
        """Retrieve PubPeer comments from API"""
        entries = []
        page = 1
        if self.last_consume_at:
            from_date = self.last_consume_at.timestamp() - 3600 * 24
        else:
            from_date = time.time() - 3600 * 24 * settings.CONS_PUBPEER_INIT_PAST
        while True:
            time.sleep(2)
            query = '{url}{page}?devkey={key}'.format(
                url=self.URL_QUERY,
                page=page,
                key=self.API_KEY
            )
            resp = requests.get(query)
            entries += json.loads(resp.text)['publications']

            ds = min([float(c['date']) for d in entries for c in d['comments']])
            if ds < from_date:
                break
            page += 1

        return entries

    def populate(self):
        """Populate DB with PubPeer comments"""
        # consume
        entries = self.consume()

        # save to database
        count = 0
        for entry in entries:
            item = self.parser.parse(entry)
            thread = self.add_or_update_entry(item)
            if thread:
                count += 1
        self.last_consume_at = datetime.datetime.now()
        return count

    def get_or_create_related_paper(self, doi):
        try:
            paper = Paper.objects.get(id_doi=doi)
        except Paper.DoesNotExist:
            paper_template = PaperParser.paper_template.copy()
            paper_template['id_doi'] = doi
            entry = {'paper': paper_template}
            new_entry = Consolidate(entry).consolidate()
            new_entry['paper']['source'] = 'PPR'
            new_entry['is_trusted'] = True
            paper, _ = PaperManager().get_or_create_from_entry(new_entry)
        return paper

    def add_or_update_entry(self, item):

        thread_entry = item['thread']
        pubpeer_entry = item['pubpeer']
        comments_entry = item['comments']

        # Thread
        try:
            thread = Thread.objects.get(pubpeer__pubpeer_id=
                                        pubpeer_entry['pubpeer_id'])
        except Thread.DoesNotExist:
            doi = pubpeer_entry['doi']
            paper = self.get_or_create_related_paper(doi)
            thread_entry['paper'] = paper.id
            thread_entry['title'] = 'Comment on: {0}'.format(paper.title)
            thread_entry['user'] = User.objects.get(
                email=settings.CONSUMER_PUBPEER_USER_EMAIL
            ).id
            form = ThreadForm(thread_entry)
            if form.is_valid():
                thread = form.save()
            else:
                raise ValueError('ThreadForm is invalid {0}'
                                 .format(form._errors))

        # PubPeer
        if thread:
            try:
                pb = PubPeer.objects.get(
                    pubpeer_id=pubpeer_entry['pubpeer_id']
                )
            except PubPeer.DoesNotExist:
                pubpeer_entry['thread'] = thread.id
                form = PubPeerForm(pubpeer_entry)
                if form.is_valid():
                    pb = form.save()
                else:
                    raise ValueError('PubPeerForm is invalid {0}'
                                     .format(form._errors))

            # PubPeerComment
            pbcs = []
            for c in comments_entry:
                try:
                    pbc = PubPeerComment.objects.get(
                        pubpeercomment_id=c['pubpeercomment_id']
                    )
                except PubPeerComment.DoesNotExist:
                    c['pubpeer'] = pb.id
                    form = PubPeerCommentForm(c)
                    if form.is_valid():
                        pbc = form.save()
                    else:
                        raise ValueError('PubPeerComment form is invalid {0}'
                                         .format(form._errors))
                pbcs.append(pbc)
