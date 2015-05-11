from django.core.exceptions import ValidationError
from stdnum import issn as issn_checker


def validate_issn(issn):
    """ Raise an Exception if issn not valid

    """
    issn_checker.validate(issn)


def validate_author_names(name):
    # process first name
    initials = [name[0] for name in name.split(' ')]
    if any(map(str.islower, initials)):
        msg = u"Name(s) must be capitalized"
        raise ValidationError(msg)
