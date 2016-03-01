# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from .constants import LANGUAGES, LANGUAGES_DETECT


def langcode_to_langpap(lang_code):
    """Convert the language code used in language detection in etalia
    language code used in the app
    """
    return dict(LANGUAGES_DETECT)[lang_code.upper()]