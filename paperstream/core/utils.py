import os
from django.core.exceptions import ImproperlyConfigured


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