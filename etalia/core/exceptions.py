# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import


class ConsolidateManagerException(Exception):
    """Top-level error for ConsolidateManager exception.
    This exception should normally not be raised, only subclasses of this
    exception."""

    def __str__(self):
        return getattr(self, 'message', '')


class PaperManagerException(Exception):
    """Top-level error for PaperManager exception.
    This exception should normally not be raised, only subclasses of this
    exception."""

    def __str__(self):
        return getattr(self, 'message', '')


class PubPeerManagerException(Exception):
    """Top-level error for PubPeerManager exception.
    This exception should normally not be raised, only subclasses of this
    exception."""

    def __str__(self):
        return getattr(self, 'message', '')


