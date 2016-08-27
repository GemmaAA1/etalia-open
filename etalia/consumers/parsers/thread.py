# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import datetime
from django.utils import timezone
from etalia.threads.models import PubPeerComment, PubPeer
from etalia.threads.forms import PubPeerForm, PubPeerCommentForm
from etalia.threads.constant import THREAD_PAPER, THREAD_PUBLIC


class PubPeerThreadParser(object):

    source = 'PPR'

    pubpeer_template = dict([(field, PubPeer._meta.get_field(field).default)
                             for field in PubPeerForm.Meta.fields])

    pubpeercomment_template = \
        dict([(field, PubPeerComment._meta.get_field(field).default)
              for field in PubPeerCommentForm.Meta.fields])

    def parse(self, entry):

        return {
            'thread': self.parse_thread(entry),
            'pubpeer': self.parse_pubpeer(entry),
            'comments': self.parse_comments(entry)
        }

    def parse_thread(self, entry):
        thread = {}
        thread['type'] = THREAD_PAPER
        thread['privacy'] = THREAD_PUBLIC
        thread['published_at'] = self.get_datetime_of_first_comment(entry)

        return thread

    def get_datetime_of_first_comment(self, entry):
        datetimes = []
        for c in entry['comments']:
            datetimes.append(datetime.datetime.fromtimestamp(int(c['date'])))
        return timezone.make_aware(min(datetimes))

    def parse_pubpeer(self, entry):

        pubpeer = self.pubpeer_template.copy()

        pubpeer['doi'] = entry.get('doi', '')
        pubpeer['link'] = entry.get('link', '')
        pubpeer['pubpeer_id'] = entry.get('pubpeer_id')

        return pubpeer

    def parse_comments(self, entry):

        comments = []

        for c in entry['comments']:

            pubpeercomment = self.pubpeercomment_template.copy()
            pubpeercomment['body'] = c['body']
            pubpeercomment['date'] = datetime.datetime.fromtimestamp(
                int(c['date'])
            )
            pubpeercomment['pubpeercomment_id'] = c.get('id', None)
            pubpeercomment['permalink'] = c.get('permalink', '')
            pubpeercomment['rating'] = c.get('rating', 0)
            pubpeercomment['user'] = c.get('user', '')
            comments.append(pubpeercomment)

        return comments