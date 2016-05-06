# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from ..models import Paper, Journal, Author
from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin


class JournalSerializer(serializers.HyperlinkedModelSerializer):
    """Journal serializer"""

    class Meta:
        model = Journal
        extra_kwargs = {
            'link': {'view_name': 'api:journal-detail'}
        }
        fields = (
            'id',
            'link',
            'title',
            'short_title',
            'url',
        )


class AuthorSerializer(serializers.HyperlinkedModelSerializer):
    """Author serializer"""

    class Meta:
        model = Author
        extra_kwargs = {
            'link': {'view_name': 'api:author-detail'}
        }
        fields = (
            'id',
            'link',
            'first_name',
            'last_name',
        )


class PaperSerializer(One2OneNestedLinkSwitchMixin,
                      serializers.HyperlinkedModelSerializer):
    """Paper serializer"""

    class Meta:
        model = Paper
        extra_kwargs = {
            'link': {'view_name': 'api:paper-detail'},
            'journal': {'view_name': 'api:journal-detail'},
            'authors': {'view_name': 'api:author-detail'},
        }
        fields = (
            'id',
            'link',
            'id_doi',
            'id_pmi',
            'id_arx',
            'id_pii',
            'id_oth',
            'title',
            'journal',
            'authors',
            'url'
        )
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'journal': {'serializer': JournalSerializer},
        }


class PaperNestedSerializer(PaperSerializer):
    """Paper nested serializer"""

    authors = AuthorSerializer(many=True, read_only=True)

    class Meta(PaperSerializer.Meta):
        read_only_fields = (
            '__all__',
        )
