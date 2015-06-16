from django import forms
from django.contrib.auth import get_user_model
from .models import Affiliation
from .validators import validate_first_name, validate_last_name

User = get_user_model()

class UserBasicForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserBasicForm, self).__init__(*args, **kwargs)
        self.fields['first_name'].validators.append(validate_first_name)
        self.fields['last_name'].validators.append(validate_last_name)

    def clean_first_name(self):
        first_name = self.cleaned_data['first_name']
        return first_name.strip()

    def clean_last_name(self):
        last_name = self.cleaned_data['last_name']
        return last_name.strip()

    class Meta:
        model = User
        fields = ('first_name', 'last_name', 'email')

class UserBasicNoEmailForm(UserBasicForm):

    class Meta:
        model = User
        fields = ('first_name', 'last_name')


class UserAffiliationForm(forms.ModelForm):

    class Meta:
        model = Affiliation
        fields = ('department', 'institution', 'city', 'state', 'country')


class UserProfileForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('first_name', 'last_name', 'email', 'affiliation')