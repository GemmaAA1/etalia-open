import random
from stdnum import issn


def gen_issn():
    digit7 = random.randint(0, 1e7-1)
    issn_ = '{0:07d}{1}'.format(
        digit7,
        issn.calc_check_digit('{0:07d}'.format(digit7)))
    return issn_

