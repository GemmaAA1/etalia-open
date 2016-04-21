# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db import models

from django.contrib.postgres.fields import ArrayField
from django.conf import settings


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


# ------------------------------------------------------------------------------
# THREADS/NLP RELATED
# ------------------------------------------------------------------------------
# Models
# ------------------------------------------------------------------------------
class ThreadVectors(TimeStampedModel):
    """ Table Thread - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use setter/getter functions to set and get vector
    """

    thread = models.ForeignKey('etalia.threads.models.Thread', related_name='vectors')

    model = models.ForeignKey('etalia.nlp.models.Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('thread', 'model')

    def __str__(self):
        return '{thread_pk}/{name}'.format(paper_pk=self.thread.pk,
                                           name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]
