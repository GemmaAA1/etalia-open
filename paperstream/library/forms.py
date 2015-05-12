from django import forms
from .models import Journal, Paper, Author, Publisher
from .validators import validate_issn

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

    def clean_id_eissn(self):
        eissn = self.cleaned_data['id_eissn']
        validate_issn(eissn)

    def clean_id_issn(self):
        issn = self.cleaned_data['id_issn']
        validate_issn(issn)
        return issn

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


# class PaperForm(forms.ModelForm):
#
#     class Meta:
#         model = Paper
#         fields = ()