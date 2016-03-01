# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class MyRouter(object):

    def route_for_task(self, task, args=None, kwargs=None):
        if task.startswith('etalia.consumers'):
            return {'queue': 'consumers',
                    'routing_key': 'consumers'}
        if task.startswith('etalia.users'):
            return {'queue': 'default',
                    'routing_key': 'default.users'}
        if task.startswith('etalia.nlp.tasks.mostsimilar'):
            return {'queue': 'mostsimilar',
                    'routing_key': 'mostsimilar'}
        if task.startswith('etalia.feeds'):
            return {'queue': 'feed',
                    'routing_key': 'feed'}
        if task.startswith('etalia.nlp'):
            return {'queue': 'nlp',
                    'routing_key': 'nlp'}
        if task.startswith('etalia.altmetric'):
            return {'queue': 'altmetric',
                    'routing_key': 'altmetric'}
        return None
