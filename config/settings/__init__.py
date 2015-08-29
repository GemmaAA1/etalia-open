# -*- coding: utf-8 -*-

import os
from paperstream.core.utils import makedirs_p

# creating logs directory if not exist
if not(os.path.isdir('logs')):
    makedirs_p('logs')
