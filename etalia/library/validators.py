# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.core.exceptions import ValidationError
from stdnum import issn as issn_checker
from stdnum.exceptions import InvalidChecksum, InvalidFormat, InvalidLength

def validate_issn(issn):
    """ Raise an ValidationException if issn not valid
    """
    if issn:
        try:
            issn_checker.validate(issn)
        except InvalidChecksum:
            raise ValidationError(
                'ISSN invalid checksum: %(issn)s',
                code='invalid',
                params={'issn': issn})
        except InvalidFormat:
            raise ValidationError(
                'ISSN invalid format: %(issn)s',
                code='invalid',
                params={'issn': issn})
        except InvalidLength:
            raise ValidationError(
                'ISSN invalid length: %(issn)s',
                code='invalid',
                params={'issn': issn})


def validate_author_names(name):
    # process first name
    initials = [name[0] for name in name.split(' ')]
    if any(map(str.islower, initials)):
        msg = u"Name(s) must be capitalized"
        raise ValidationError(msg)
