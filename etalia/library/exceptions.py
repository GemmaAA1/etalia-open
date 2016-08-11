# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.utils.encoding import force_text


class PubmedException(Exception):
    """
    Base class for PubMed exceptions.
    """
    default_detail = 'A Pubmed error occurred.'

    def __init__(self, detail=None):
        if detail is not None:
            self.detail = force_text(detail)
        else:
            self.detail = force_text(self.default_detail)

    def __str__(self):
        return self.detail


class CrossrefException(Exception):
    """
    Base class for Crossref exceptions.
    """
    default_detail = 'A CrossRef error occurred.'

    def __init__(self, detail=None):
        if detail is not None:
            self.detail = force_text(detail)
        else:
            self.detail = force_text(self.default_detail)

    def __str__(self):
        return self.detail


class ArxivException(Exception):
    """
    Base class for Arxiv exceptions.
    """
    default_detail = 'An Arxiv error occurred.'

    def __init__(self, detail=None):
        if detail is not None:
            self.detail = force_text(detail)
        else:
            self.detail = force_text(self.default_detail)

    def __str__(self):
        return self.detail