# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


# List of tuples that defines task routing
# Task name is compared to first argument of tuple with 'startswith'. First
# occurrence that match is used for routing


class MyRouter(object):

    TASKS_ROUTING_MAP = [
        ('etalia.consumers', {'queue': 'consumers', 'routing_key': 'consumers'}),
        ('etalia.users', {'queue': 'default', 'routing_key': 'default.users'}),
        ('etalia.nlp.tasks.pe_dispatcher', {'queue': 'pe', 'routing_key': 'pe'}),
        ('etalia.nlp.tasks.te_dispatcher', {'queue': 'te', 'routing_key': 'te'}),
        ('etalia.nlp.tasks.nlp_dispatcher', {'queue': 'nlp', 'routing_key': 'nlp'}),
        ('etalia.nlp.tasks.paperengine', {'queue': 'update_engines', 'routing_key': 'update.pe'}),
        ('etalia.nlp.tasks.threadengine', {'queue': 'update_engines', 'routing_key': 'update.te'}),
        ('etalia.nlp', {'queue': 'nlp', 'routing_key': 'nlp'}),
        ('etalia.feeds', {'queue': 'feed', 'routing_key': 'feed'}),
        ('etalia.altmetric', {'queue': 'altmetric', 'routing_key': 'altmetric'}),
    ]

    def route_for_task(self, task, args=None, kwargs=None):
        for pattern, route in self.TASKS_ROUTING_MAP:
            if task.startswith(pattern):
                return route.copy()
        return None