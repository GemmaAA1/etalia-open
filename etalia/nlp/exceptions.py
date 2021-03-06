# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class ModelError(Exception):
    """Top-level error for model exception.
    This exception should normally not be raised, only subclasses of this
    exception."""

    def __str__(self):
        return getattr(self, 'message', '')


class InvalidState(ModelError):
    """The state of the model is wrong
    This generally means that the object is currently doing something else
    and/or the current task cannot be processed.
    """

    message = 'Current state is invalid'

    def __init__(self, message):
        super(InvalidState, self).__init__(message)
        self.message = message
