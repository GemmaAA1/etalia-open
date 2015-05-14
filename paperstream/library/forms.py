from django import forms
from .models import Journal, Paper, Author, Publisher, CorpAuthor
import requests
from requests.exceptions import RequestException


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


class JournalFormFillUp(JournalForm):
    """This form is use to add new data from API, not from view
    (at least originally, e.g. in populate app)
    The Form has the following behavior:
        - if the form is not instantiated with a Journal instance, it behaves
        as JournalForm
        - if the form is instantiated with a Journal instance:
            -- if new data value is blank, it keeps the initial value
            -- if new data value is not blank, it replaced initial value
    """

    def __init__(self, *args, **kwargs):
        super(JournalFormFillUp, self).__init__(*args, **kwargs)

    # clean fields with instance fields
    def clean_id_issn(self):
        return self.fill_up('id_issn')

    def clean_id_arx(self):
        return self.fill_up('id_arx')

    def clean_id_eissn(self):
        return self.fill_up('id_eissn')

    def clean_id_oth(self):
        return self.fill_up('id_oth')

    def clean_title(self):
        return self.fill_up('title')

    def clean_short_title(self):
        return self.fill_up('short_title')

    def clean_period(self):
        return self.fill_up('period')

    def clean_url(self):
        return self.fill_up('url')

    def clean_scope(self):
        return self.fill_up('scope')

    def clean_publisher(self):
        return self.fill_up('publisher')

    def fill_up(self, field):
        # call super clean_<field> method
        super_method = getattr(super(JournalFormFillUp, self),
                               'clean_' + field, None)

        if super_method:
            self.cleaned_data[field] = super_method()

        field_val = self.cleaned_data[field]
        # fill up blank or erase initial value if not blank
        if field_val:
            return field_val
        elif self.instance.pk:
            return getattr(self.instance, field)
        else:
            return field_val


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
        fields = ('type',
                  'id_doi',
                  'id_arx',
                  'id_pii',
                  'id_pmi',
                  'id_oth',
                  'title',
                  'abstract',
                  'journal',
                  'volume',
                  'issue',
                  'page',
                  'date_ep',
                  'date_p',
                  'date_lr',
                  'url',
                  'language',
                  'is_aip',
                  'is_pre_print',
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
                self.Meta.model.detects_language(text)

        return cleaned_data
