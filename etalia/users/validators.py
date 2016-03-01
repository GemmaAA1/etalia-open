# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from __future__ import unicode_literals
from etalia.core.validators import CustomRegexValidator

def validate_first_name(char):
    regex_validator = CustomRegexValidator(
        regex=r'^[\p{L}\s\']+$',
        message="First name must contain only letters",
        code='invalid')
    regex_validator(char)

def validate_last_name(char):
    regex_validator = CustomRegexValidator(
        regex=r'^[\p{L}\s\']+$',
        message="Last name must contain only letters",
        code='invalid')
    regex_validator(char)

