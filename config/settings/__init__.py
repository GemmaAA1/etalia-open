# -*- coding: utf-8 -*-

import os
import re
from paperstream.core.utils import makedirs_p

# creating logs directory if not exist
if not(os.path.isdir('logs')):
    makedirs_p('logs')


def get_version(package):
    """
    Return package version as listed in `__version__` in `init.py`.
    """
    init_py = open(os.path.join(package, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)