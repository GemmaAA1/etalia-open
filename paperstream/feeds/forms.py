from django import forms
from .models import UserFeed

class CreateUserFeedForm(forms.ModelForm):

    class Meta:
        model = UserFeed
        fields = ('name', )
        widgets = {
            'name':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Name',
                                                'initial': ''}),
        }

    def __init__(self, *args, **kwargs):
        self.request = None
        if 'request' in kwargs:
            self.request = kwargs.pop('request')
        super(CreateUserFeedForm, self).__init__(*args, **kwargs)
        self.fields['name'].initial = ''

    def clean(self):
        cleaned_data = self.cleaned_data
        if self.request:        # add user to form data
            cleaned_data['user'] = self.request.user
        return cleaned_data


class CheckedForm(forms.Form):

    checked = forms.BooleanField()

