# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class ReferenceManagerInstanceError(Error):
    """Exception raised for error on social_auth reference manager.

    Args:
        user (User): User instance
        message (str): explanation of the error
    """

    def __init__(self, user, message):
        self.user = user
        self.message = message


class OrcidInstanceError(Error):
    """Exception raised for error on social_auth Orcid.

    Args:
        user (User): User instance
        message (str): explanation of the error
    """

    def __init__(self, user, message):
        self.user = user
        self.message = message