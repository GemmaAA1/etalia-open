import logging

logger = logging.getLogger(__name__)


def test():
    for a in range(10):
        logger.info('test %d' % a)