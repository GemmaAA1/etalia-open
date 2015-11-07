# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.contrib.auth import get_user_model
from django.contrib.auth.forms import AuthenticationForm
from .models import Affiliation, UserSettings
from .validators import validate_first_name, validate_last_name

from paperstream.nlp.models import Model
from paperstream.feeds.constants import STREAM_METHODS, TREND_METHODS

User = get_user_model()


class UserAuthenticationForm(AuthenticationForm):

    error_messages = {
        'invalid_login': "Please enter a correct %(username)s and password.",
        'inactive': "This account is inactive.",
    }


class UserBasicForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('email', 'first_name', 'last_name')

        widgets = {
            'email':
                forms.TextInput(attrs={'class': 'form-control input-lg ',
                                                'placeholder': 'Email'}),
            'first_name':
                forms.TextInput(attrs={'class': 'form-control input-lg ',
                                                'placeholder': 'First Name'}),
            'last_name':
                forms.TextInput(attrs={'class': 'form-control input-lg ',
                                                'placeholder': 'Last Name'})
        }

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


class UpdateUserNameForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UpdateUserNameForm, self).__init__(*args, **kwargs)
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
        fields = ('first_name', 'last_name')
        widgets = {
            'first_name':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'First Name'}),
            'last_name':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Last Name'}),
        }


class UpdateUserTitleForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('title', )
        widgets = {
            'title':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Title'}),
        }


class UpdateUserPositionForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('position', )
        widgets = {
            'position':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Position'}),
        }


class UserAffiliationForm(forms.ModelForm):

    class Meta:
        model = Affiliation
        fields = ('department', 'institution', 'city', 'state', 'country')
        widgets = {
            'department':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Department'}),
            'institution':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Institution'}),
            'city':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'City'}),
            'state':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'State'}),
            'country':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Country'}),
        }


class UserStreamSettingsForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserStreamSettingsForm, self).__init__(*args, **kwargs)
        self.fields['stream_model'].choices = [(mod.pk, mod.name)
                                        for mod in Model.objects.filter(is_active=True)]

    class Meta:
        model = UserSettings
        fields = ('stream_model', 'stream_time_lapse', 'stream_method',
                  )
        widgets = {
            'stream_model': forms.Select(attrs={'class': 'form-control'}),
            'stream_time_lapse': forms.Select(attrs={'class': 'form-control'}),
            'stream_method': forms.Select(attrs={'class': 'form-control'}),
        }


class UserTrendSettingsForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserTrendSettingsForm, self).__init__(*args, **kwargs)
        self.fields['trend_model'].choices = [(mod.pk, mod.name)
                                        for mod in Model.objects.filter(is_active=True)]

    class Meta:
        model = UserSettings
        fields = ('trend_model', 'trend_time_lapse', 'trend_method'
                  )
        widgets = {
            'trend_model': forms.Select(attrs={'class': 'form-control'}),
            'trend_time_lapse': forms.Select(attrs={'class': 'form-control'}),
            'trend_method': forms.Select(attrs={'class': 'form-control'}),
        }
