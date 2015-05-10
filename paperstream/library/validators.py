from django.core.exceptions import ValidationError
from stdnum import issn as issn_checker


def validate_id_issn(issn):
    """ Raise an Exception if issn not valid

    """
    issn_checker.validate(issn)

def validate_id_eissn(eissn):
    """ Raise an Exception if eissn not valid
    """
    validate_id_issn(eissn)

