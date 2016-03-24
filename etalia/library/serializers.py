# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers

from .models import Paper, Journal, Author


class JournalSerializer(serializers.ModelSerializer):

    class Meta:
        model = Journal
        fields = ('id',
                  'title',
                  'short_title',
                  'url')


class AuthorSerializer(serializers.ModelSerializer):

    class Meta:
        model = Author
        fields = ('id',
                  'first_name',
                  'last_name',)


class PaperSerializer(serializers.ModelSerializer):

    authors = AuthorSerializer(many=True, read_only=True)
    journal = JournalSerializer(many=False, read_only=True)

    class Meta:
        model = Paper
        fields = ('id',
                  'title',
                  'journal',
                  'authors',
                  'id_doi',
                  'id_pmi',
                  'id_arx',
                  'url')




