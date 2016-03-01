# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import regex as reg
from django.utils import six
from django.core.validators import RegexValidator


class CustomRegexValidator(RegexValidator):
    def __init__(self, regex=None, message=None, code=None, inverse_match=None, flags=None):
        if regex is not None:
            self.regex = regex
        if message is not None:
            self.message = message
        if code is not None:
            self.code = code
        if inverse_match is not None:
            self.inverse_match = inverse_match
        if flags is not None:
            self.flags = flags
        if self.flags and not isinstance(self.regex, six.string_types):
            raise TypeError("If the flags are set, regex must be a regular expression string.")

        if isinstance(self.regex, six.string_types):
            self.regex = reg.compile(self.regex, self.flags)