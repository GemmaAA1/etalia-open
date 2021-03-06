# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import sys
import traceback
import numpy as np
from django.core import mail
from django.views.debug import ExceptionReporter
from django.core.exceptions import ValidationError
from django.conf import settings
from django.db import connection


class AttrDict(dict):
    """To access dictionary (key,value) as attribute"""
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def get_celery_worker_status():
    ERROR_KEY = "ERROR"
    try:
        from celery.task.control import inspect
        insp = inspect()
        d = insp.stats()
        if not d:
            d = { ERROR_KEY: 'No running Celery workers were found.' }
    except IOError as e:
        from errno import errorcode
        msg = "Error connecting to the backend: " + str(e)
        if len(e.args) > 0 and errorcode.get(e.args[0]) == 'ECONNREFUSED':
            msg += ' Check that the RabbitMQ server is running.'
        d = {ERROR_KEY: msg}
    except ImportError as e:
        d = {ERROR_KEY: str(e)}
    return d


def pad_or_trim_vector(vector):
    if isinstance(vector, np.ndarray):
        if vector.size == 1:
            vector = [np.squeeze(vector).tolist()]
        else:
            vector = np.squeeze(vector).tolist()

    if len(vector) <= settings.NLP_MAX_VECTOR_SIZE:
        vector += [0.0] * (settings.NLP_MAX_VECTOR_SIZE - len(vector))
    else:
        return vector[:settings.NLP_MAX_VECTOR_SIZE]

    return vector


def pad_neighbors(vector):
    if isinstance(vector, np.ndarray):
        if vector.size == 1:
            vector = [np.squeeze(vector).tolist()]
        else:
            vector = np.squeeze(vector).tolist()

    if len(vector) <= settings.NLP_MAX_KNN_NEIGHBORS:
        vector += [0] * (settings.NLP_MAX_KNN_NEIGHBORS - len(vector))
    else:
        raise ValidationError('vector is larger than NLP_MAX_KNN_NEIGHBORS')
    return vector


def makedirs_p(path):
    if sys.version_info < (3, 4):
        try:
            os.makedirs(path)
        except OSError:
            if not os.path.isdir(path):
                raise
    else:
        os.makedirs(path, exist_ok=True)


def db_table_exists(table_name):
    return table_name in connection.introspection.table_names()