# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models

# Create your models here.

# Abstract models
class TimeStampedModel(models.Model):
    """
    An abstract base class model that provides selfupdating
    ``created`` and ``modified`` fields.
    """
    created = models.DateTimeField(auto_now_add=True)
    modified = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True


# Custom fields
class NullableCharField(models.CharField):
    description = "CharField that stores NULL but returns ''"
    __metaclass__ = models.SubfieldBase

    def to_python(self, value):
        if isinstance(value, models.EmailField):
            return value
        return value or ''

    def get_prep_value(self, value):
        return value or None
