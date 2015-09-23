# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import requests
from requests.exceptions import RequestException
from django import forms

from .models import Journal, Paper, Author, Publisher, CorpAuthor
from .constants import SOURCE_TYPE


class FillBlanksMixin(object):
    """This Mixin is use to give the form the following behavior:
        - if the form is not instantiated with a Journal instance, it behaves
        as JournalForm
        - if the form is instantiated with a Journal instance:
            -- if new data value is blank, it keeps the initial value
            -- if new data value is not blank, it replaced initial value
    """

    def clean(self):
        """Order matters. fill_blank, which is field related, needs to be called
        before super().clean(). Field wise cleaning is done here for DRY.
        :return: (dict)
        """
        for field in self.cleaned_data.keys():
            self.cleaned_data[field] = self.fill_blank(field)

        cleaned_data = super(FillBlanksMixin, self).clean()
        return cleaned_data

    def fill_blank(self, field):
        """Update field but do not update with blank or None value
        :param field (str):
        :return:
        """
        field_val = self.cleaned_data[field]
        if field_val:
            return field_val
        elif self.instance.pk:
            return getattr(self.instance, field)
        else:
            return field_val


class PublisherForm(forms.ModelForm):

    class Meta:
        model = Publisher
        fields = ('name', 'url')


class JournalForm(forms.ModelForm):
    class Meta:
        model = Journal
        fields = ('id_issn',
                  'id_eissn',
                  'id_arx',
                  'id_oth',
                  'title',
                  'short_title',
                  'period',
                  'url',
                  'scope',
                  'publisher')

    def __init__(self, *args, **kwargs):
        super(JournalForm, self).__init__(*args, **kwargs)
        # publisher foreign key is defined by the publisher name (unique)
        self.fields['publisher'] = \
            forms.ModelChoiceField(queryset=Publisher.objects.all(),
                                   empty_label=None,
                                   to_field_name='name',
                                   required=False)

    def clean_title(self):
        title = self.cleaned_data['title']
        # if title all upper or lower, capitalized
        if title.isupper() or title.islower():
            title = title.title()
        return title


class JournalFormFillBlanks(FillBlanksMixin, JournalForm):
    pass


class AuthorForm(forms.ModelForm):

    class Meta:
        model = Author
        fields = ('first_name',
                  'last_name',
                  'email')

    def clean_first_name(self):
        # capitalize
        return self.cleaned_data['first_name'].title()

    def clean_last_name(self):
        # capitalize
        return self.cleaned_data['last_name'].title()


class CorpAuthorForm(forms.ModelForm):

    class Meta:
        model = CorpAuthor
        fields = ('name',
                  )

    def clean_name(self):
        # capitalize
        return self.cleaned_data['name'].title()


class PaperForm(forms.ModelForm):

    class Meta:
        model = Paper
        fields = ('publish_status',
                  'type',
                  'id_doi',
                  'id_arx',
                  'id_pii',
                  'id_pmi',
                  'id_isbn',
                  'id_oth',
                  'title',
                  'abstract',
                  'journal',
                  'volume',
                  'issue',
                  'page',
                  'date_ep',
                  'date_pp',
                  'date_lr',
                  'url',
                  'language',
                  'source',
                  )

    def __init__(self, *args, **kwargs):
        super(PaperForm, self).__init__(*args, **kwargs)
        # publisher foreign key is defined by the publisher name (unique)
        self.fields['journal'] = \
            forms.ModelChoiceField(queryset=Journal.objects.all(),
                                   empty_label=None,
                                   to_field_name='title',
                                   required=False)
        self.fields['source'] = \
            forms.ChoiceField(choices=SOURCE_TYPE)

    def clean_title(self):
        """Removing trailling/heading brackets found in some title
        """
        title = self.cleaned_data['title']
        if (title[0], title[-1]) == ('[', ']'):
            title = title[1:-1]
        return title

    def clean_abstract(self):
        abstract = self.cleaned_data['abstract']
        # remove word 'abstract' if it is first word
        if abstract.lower().startswith('abstract') \
                or abstract.lower().startswith('abstract:'):
            abstract = abstract.lower().split(' ', 1)[1]

        return abstract

    def clean(self):
        cleaned_data = super(PaperForm, self).clean()

        # language
        if not self.cleaned_data['language']:
            text = ' '.join([self.cleaned_data.get('abstract', ''),
                             self.cleaned_data.get('title', '')])
            if text.strip():
                self.cleaned_data['language'] = \
                    self.Meta.model.detect_language(text)

        return cleaned_data


class PaperFormFillBlanks(FillBlanksMixin, PaperForm):
    pass


class PaperFormClean(FillBlanksMixin, PaperForm):

    def clean_url(self):
        # Check if URL returns 200
        url = self.cleaned_data['url']
        if url:
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    return url
                else:
                    return ''
            except RequestException:
                return ''

    def clean_id_doi(self):
        """Check if doi exist of doi.org"""
        # format
        id_doi = self.cleaned_data['id_doi'].lower()

        # Check if doi valid requesting http://doi.org/<doi>
        if id_doi:
            url = 'http://doi.org/{doi}'.format(doi=id_doi)
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    return id_doi
                else:
                    return ''
            except RequestException:
                return ''

    # def clean_id_arx(self):
    #     #TODO:
    #     raise NotImplemented
    #
    # def clean_id_pii(self):
    #     #TODO:
    #     raise NotImplemented
    #
    # def clean_id_pmi(self):
    #     #TODO:
    #     raise NotImplemented
    #
    # def clean_journal(self):
    #     #TODO:
    #     raise NotImplemented

    def clean(self):
        cleaned_data = super(PaperForm, self).clean()

        # Detect language
        if not self.cleaned_data['language']:
            text = ' '.join([self.cleaned_data['abstract'],
                             self.cleaned_data['title']])
            self.cleaned_data['language'] = \
                self.Meta.model.detect_language(text)

        # Test URL
        # try with id_doi
        if not cleaned_data['url']:
            id_doi = cleaned_data.get('id_doi', '')
            if id_doi:
                url = 'http://doi.org/{doi}'.format(doi=id_doi)
                try:
                    response = requests.get(url)
                    if response.status_code == 200:
                        cleaned_data['url'] = url
                except RequestException:
                    pass
        # try with id_pmi
        if not cleaned_data['url']:
            id_pmi = cleaned_data.get('id_pmi', '')
            if id_pmi:
                url = 'http://www.ncbi.nlm.nih.gov/pubmed/{pmi}'.format(
                    pmi=id_pmi)
                try:
                    response = requests.get(url)
                    if response.status_code == 200:
                        cleaned_data['url'] = url
                except RequestException:
                    pass

        return cleaned_data
