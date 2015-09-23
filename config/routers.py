# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class MyRouter(object):

    def route_for_task(self, task, args=None, kwargs=None):
        if task.startswith('paperstream.consumers'):
            return {'queue': 'consumers',
                    'routing_key': 'consumers'}
        if task.startswith('paperstream.users'):
            return {'queue': 'default',
                    'routing_key': 'default.users'}
        if task == 'paperstream.nlp.models.populate_neighbors':
            return {'queue': 'lsh',
                    'routing_key': 'lsh.populate_neighbors'}
        if task.startswith('paperstream.feeds'):
            return {'queue': 'lsh',
                    'routing_key': 'lsh.feeds'}
        if task.startswith('paperstream.nlp'):
            return {'queue': 'nlp',
                    'routing_key': 'nlp'}
        if task.startswith('paperstream.altmetric'):
            return {'queue': 'altmetric',
                    'routing_key': 'altmetric'}
        return None
