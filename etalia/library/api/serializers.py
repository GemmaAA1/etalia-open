# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from rest_framework.reverse import reverse

from ..models import Paper, Journal, Author, PaperUser
from ..constants import PAPER_ADDED, PAPER_TRASHED, PAPER_STORE
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

    state = serializers.SerializerMethodField()

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
            'abstract',
            'url',
            'state',
        )
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'journal': {'serializer': JournalSerializer},
        }

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        paperuser = obj.state(self.context['request'].user)
        if paperuser:
            if self.one2one_nested:
                return PaperUserSerializer(
                    instance=paperuser,
                    context={'request': self.context['request']},
                    one2one_nested=False
                ).data
            else:
                return reverse('api:paperuser-detail',
                               kwargs={'pk': paperuser.id},
                               request=self.context.get('request', None),
                               format=self.context.get('format', None))
        return None


class PaperNestedSerializer(PaperSerializer):
    """Paper nested serializer"""

    authors = AuthorSerializer(many=True, read_only=True)

    class Meta(PaperSerializer.Meta):
        read_only_fields = (
            '__all__',
        )


class JournalFilterSerializer(serializers.ModelSerializer):
    """Journal serializer in filter side panel"""

    count = serializers.IntegerField(read_only=True)
    label = serializers.SerializerMethodField()

    class Meta:
        model = Journal
        fields = (
            'id',
            'label',
            'count')
        read_only_fields = (
            '__all__'
        )

    def get_label(self, obj):
        return '{0}'.format(obj.title)


class AuthorFilterSerializer(serializers.ModelSerializer):
    """Journal serializer in filter side panel"""

    count = serializers.IntegerField(read_only=True)
    label = serializers.SerializerMethodField()

    class Meta:
        model = Author
        fields = (
            'id',
            'label',
            'count')
        read_only_fields = (
            '__all__'
        )

    def get_label(self, obj):
        return '{0} {1}'.format(obj.first_name, obj.last_name)


class PaperFilterSerializer(serializers.BaseSerializer):
    """Serializer for filters on side panel of Threads list"""

    def to_representation(self, instance):
        return {
            'groups': [
                {
                    "name": "author_id",
                    "label": "Authors",
                    "entries": [
                        AuthorFilterSerializer(instance=author,
                                               context=self.context).data
                        for author in instance.get('authors', None)
                        ]
                },
                {
                    "name": "journal_id",
                    "label": "Journals",
                    "entries": [
                        JournalFilterSerializer(instance=journal,
                                                context=self.context).data
                        for journal in instance.get('journals', None)
                        ]
                }
            ]
        }


class PaperUserSerializer(One2OneNestedLinkSwitchMixin,
                          serializers.HyperlinkedModelSerializer):

    class Meta:
        model = PaperUser
        extra_kwargs = {
            'link': {'view_name': 'api:paperuser-detail'},
            'user': {'view_name': 'api:user-detail'},
            'paper': {'view_name': 'api:paper-detail'},
        }
        fields = (
            'id',
            'link',
            'user',
            'paper',
            'watch',
            'store'
        )
        read_only_fields = (
            'id',
            'link',
        )
        switch_kwargs = {
            # 'user': {'serializer': UserSerializer,
            #          'one2one_nested': False},
            'paper': {'serializer': PaperSerializer,
                      'one2one_nested': False},
        }

    def validate_user(self, value):
        """User passed is same as user login"""
        if not value == self.context['request'].user:
            raise serializers.ValidationError(
                "user deserialized and user logged are different")
        return value

    def validate_store(self, value):
        """Check that paper is added if request is to trashed paper"""
        if self.instance:
            if value == PAPER_TRASHED and \
                    not self.instance.store == PAPER_ADDED:
                raise serializers.ValidationError(
                    "You cannot trash a paper that you have not added")
        return value


class PaperUserUpdateSerializer(PaperUserSerializer):

    class Meta(PaperUserSerializer.Meta):

        read_only_fields = (
            'id',
            'link',
            'paper',
            'user'
        )

    def update(self, instance, validated_data):
        """Trigger add() and trash() method as needed"""
        serializers.raise_errors_on_nested_writes('update', self, validated_data)
        err = None
        for attr, value in validated_data.items():
            if attr == 'store':
                if value == PAPER_ADDED:
                    err = instance.add()
                elif value == PAPER_TRASHED:
                    err = instance.trash()
                else:
                    raise serializers.ValidationError(
                        'store value outside of choices ({0})'.format(PAPER_STORE))
                if err:
                    raise serializers.ValidationError(err)
            else:
                setattr(instance, attr, value)
        instance.save()

        return instance
