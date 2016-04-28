# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from ..models import Paper, Journal, Author


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


class PaperSerializer(serializers.HyperlinkedModelSerializer):
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


class PaperNestedSerializer(serializers.HyperlinkedModelSerializer):
    """Paper nested serializer"""

    authors = AuthorSerializer(many=True, read_only=True)
    journal = JournalSerializer(many=False, read_only=True)

    class Meta:
        model = Paper
        extra_kwargs = {
            'link': {'view_name': 'api:paper-detail'},
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
