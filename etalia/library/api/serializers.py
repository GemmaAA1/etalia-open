# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from mendeley.exception import MendeleyApiException

from rest_framework import serializers
from rest_framework.reverse import reverse

from django.utils.text import slugify
from etalia.core.api.mixins import One2OneNestedLinkSwitchMixin
from etalia.core.api.exceptions import MendeleyRedirectLoginErrorSerializer
from etalia.feeds.models import StreamPapers, TrendPapers
from etalia.users.models import UserLibPaper


from ..models import Paper, Journal, Author, PaperUser
from ..constants import PAPER_ADDED, PAPER_TRASHED, PAPER_STORE
from ..mixins import PaperEagerLoadingMixin


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


class PaperSerializer(PaperEagerLoadingMixin,
                      One2OneNestedLinkSwitchMixin,
                      serializers.HyperlinkedModelSerializer):
    """Paper serializer"""

    state = serializers.SerializerMethodField()
    linked_threads_count = serializers.SerializerMethodField()
    authors = serializers.SerializerMethodField()
    url = serializers.SerializerMethodField()
    slug = serializers.SerializerMethodField()
    new = serializers.SerializerMethodField()
    is_orcid = serializers.SerializerMethodField()

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
            # 'date',
            'date_co',
            'date_ep',
            'date_pp',
            'slug',
            'linked_threads_count',
            'state',
            'new',
            'is_orcid',
        )
        read_only_fields = (
            '__all__',
        )
        switch_kwargs = {
            'journal': {'serializer': JournalSerializer},
        }

    def get_slug(self, obj):
        return '{slug}_{id}'.format(slug=slugify(obj.title),
                                    id=obj.id)

    def get_url(self, obj):
        if obj.id_doi:
            return obj.get_doi_url()
        else:
            return obj.url

    def get_authors(self, obj):
        id_aut = dict([(a.id, a) for a in obj.authors.all()])
        pos_id = [(ap.position, ap.author_id) for ap in obj.authorpaper_set.all()]
        pos_id = sorted(pos_id, key=lambda x: x[0])
        authors = [id_aut[x[1]] for x in pos_id]
        return [reverse('api:author-detail',
                        kwargs={'pk': auth.id},
                        request=self.context.get('request', None),
                        format=self.context.get('format', None))
                for auth in authors]

    def get_state(self, obj):
        """Get state based on ThreadUser instance if exists"""
        if hasattr(obj, 'pu') and obj.pu:
            if self.one2one_nested:
                return PaperUserSerializer(
                    instance=obj.pu[0],
                    context={'request': self.context['request']},
                    one2one_nested=False
                ).data
            else:
                return reverse('api:paperuser-detail',
                               kwargs={'pk': obj.pu.id},
                               request=self.context.get('request', None),
                               format=self.context.get('format', None))
        return None

    def get_new(self, obj):
        request = self.context['request']
        scored = request.query_params.get('scored')
        type = request.query_params.get('type', 'stream')
        feed_name = request.query_params.get('feed', 'main')
        try:
            if scored == '1':
                if type == 'stream':
                    if hasattr(obj, 'sp') and obj.sp:
                        return obj.sp[0].new
                    else:
                        return obj.streampapers_set\
                            .get(stream__user=request.user,
                                 stream__name=feed_name)\
                            .new
                if type == 'trend':
                    if hasattr(obj, 'tp') and obj.tp:
                        return obj.tp[0].new
                    else:
                        return obj.trendpapers_set\
                            .get(trend__user=request.user,
                                 trend__name=feed_name)\
                            .new
        except StreamPapers.DoesNotExist or TrendPapers.DoesNotExist:
            pass

        return None

    def get_is_orcid(self, obj):
        request = self.context['request']
        if request.user.is_anonymous():
            return False
        else:
            try:
                return request.user.lib.userlib_paper.get(paper=obj.id).is_orcid
            except UserLibPaper.DoesNotExist:
                return False

    def get_linked_threads_count(self, obj):
        if hasattr(obj, 'threads'):
            return len(obj.threads)
        else:
            return None


class PaperNestedSerializer(PaperSerializer):
    """Paper nested serializer"""

    authors = serializers.SerializerMethodField()

    class Meta(PaperSerializer.Meta):
        read_only_fields = (
            '__all__',
        )

    def get_authors(self, obj):
        id_aut = dict([(a.id, a) for a in obj.authors.all()])
        pos_id = [(ap.position, ap.author_id) for ap in obj.authorpaper_set.all()]
        pos_id = sorted(pos_id, key=lambda x: x[0])
        authors = [id_aut[x[1]] for x in pos_id]
        return AuthorSerializer(
            authors,
            many=True,
            context={'request': self.context.get('request', None)})\
            .data


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
        if obj.short_title and (len(obj.short_title) > 1 or len(obj.title) == 1):
            return '{0}'.format(obj.short_title)
        else:
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
                    "name": "journal_id",
                    "label": "Journals",
                    "entries": [
                        JournalFilterSerializer(instance=journal,
                                                context=self.context).data
                        for journal in instance.get('journals', None)
                        ]
                },
                {
                    "name": "author_id",
                    "label": "Authors",
                    "entries": [
                        AuthorFilterSerializer(instance=author,
                                               context=self.context).data
                        for author in instance.get('authors', None)
                        ]
                }
            ]
        }


class PaperUserSerializer(One2OneNestedLinkSwitchMixin,
                          serializers.HyperlinkedModelSerializer):

    is_orcid = serializers.SerializerMethodField()

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
            'store',
            'is_orcid'

        )
        read_only_fields = (
            'id',
            'link',
        )
        switch_kwargs = {
            # 'user': {'serializer': UserSerializer,
            #          'one2one_nested': False},
            # 'paper': {'serializer': PaperSerializer,
            #           'one2one_nested': False},
        }

    def get_is_orcid(self, obj):
        try:
            return obj.user.lib.userlib_paper.get(paper_id=obj.paper.id).is_orcid
        except UserLibPaper.DoesNotExist:
            return False

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

    def update(self, instance, validated_data):
        """Subclassing for triggering add() and trash() methods as needed"""
        try:
            serializers.raise_errors_on_nested_writes('update', self, validated_data)
            for attr, value in validated_data.items():
                if attr == 'store' and value is not None:
                    if value == PAPER_ADDED and not instance.store == PAPER_ADDED:
                        err = instance.add()
                    elif value == PAPER_TRASHED and not instance.store == PAPER_TRASHED:
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
        except MendeleyApiException:
            raise MendeleyRedirectLoginErrorSerializer()

    def create(self, validated_data):
        """Subclassing for triggering add() and trash() methods as needed"""
        try:
            instance = super(PaperUserSerializer, self).create(validated_data)
            for attr, value in validated_data.items():
                if attr == 'store' and value is not None:
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
            return instance
        except MendeleyApiException:
            raise MendeleyRedirectLoginErrorSerializer()


class JsonLdSerializer(serializers.BaseSerializer):

    def to_representation(self, instance):
        return {
            "@context": "http://schema.org",
            "@graph": [
                {
                "@id": "#issue",
                "@type": "PublicationIssue",
                "issueNumber": instance.issue or "",
                "datePublished": instance.date or "",
                "isPartOf": {
                    "@id": "#periodical",
                    "@type": [
                        "PublicationVolume",
                        "Periodical"
                    ],
                    "name": instance.journal.title if instance.journal else "",
                    "issn": [i for i in [instance.journal.id_issn, instance.journal.id_eissn]
                             if instance.journal and i],
                    "volumeNumber": instance.volume,
                    "publisher": instance.print_publisher_name
                }
                },
                {
                "@type": "ScholarlyArticle",
                "isPartOf": "#issue",
                "description": instance.abstract,
                "sameAs": instance.url,
                "pageStart": instance.page.split('-')[0],
                "pageEnd": instance.page.split('-')[1] if len(instance.page.split('-')) > 1 else '',
                "name": instance.title,
                "author": instance.print_json_ld_authors
                }
            ]
        }
