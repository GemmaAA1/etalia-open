import os
import numpy as np
from django.core.exceptions import ImproperlyConfigured
from django.core.exceptions import ValidationError
from django.conf import settings


def get_env_variable(var_name, default=None):
    """
    Get the environment variable
    """
    try:
        return os.environ[var_name]
    except KeyError:
        if default:
            return default
            pass
        else:
            error_msg = "Set the {} environment variable".format(var_name)
            raise ImproperlyConfigured(error_msg)

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

def pad_vector(vector):
    if isinstance(vector, np.ndarray):
        if vector.size == 1:
            vector = [np.squeeze(vector).tolist()]
        else:
            vector = np.squeeze(vector).tolist()

    if len(vector) <= settings.NLP_MAX_VECTOR_SIZE:
        vector += [None] * (settings.NLP_MAX_VECTOR_SIZE - len(vector))
    else:
        raise ValidationError('vector is larger than NLP_MAX_VECTOR_SIZE')

    return vector

def pad_neighbors(vector):
    if isinstance(vector, np.ndarray):
        if vector.size == 1:
            vector = [np.squeeze(vector).tolist()]
        else:
            vector = np.squeeze(vector).tolist()

    if len(vector) <= settings.NLP_MAX_KNN_NEIGHBORS:
        vector += [None] * (settings.NLP_MAX_KNN_NEIGHBORS - len(vector))
    else:
        raise ValidationError('vector is larger than NLP_MAX_KNN_NEIGHBORS')
    return vector
