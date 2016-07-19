# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.contrib.auth import get_user_model
from django.contrib.auth.forms import AuthenticationForm
from django.utils import timezone
from .models import Affiliation, UserSettings, UserLibPaper
from .validators import validate_first_name, validate_last_name
from etalia.library.models import PaperUser
from etalia.library.constants import PAPER_PINNED, PAPER_BANNED

from etalia.nlp.models import Model
from etalia.feeds.constants import STREAM_METHODS, TREND_METHODS

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


class UserFingerprintSettingsForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserFingerprintSettingsForm, self).__init__(*args, **kwargs)
        if 'instance' in kwargs:
            self.fields['fingerprint_roll_back_deltatime'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].fingerprint_roll_back_deltatime)
            delta_month = int((timezone.now().date() - kwargs['instance'].user.lib.d_oldest).days / 30) + 1
            self.fields['fingerprint_roll_back_deltatime'].widget.attrs['data-slider-max'] = delta_month

    class Meta:
        model = UserSettings
        fields = (
                  'fingerprint_roll_back_deltatime',
                  )
        widgets = {
            'fingerprint_roll_back_deltatime': forms.TextInput(attrs={
                'data-slider-id': 'fingerprint_roll_back_deltatime_slider',
                'type': 'text',
                'data-slider-min': '1',
                'data-slider-max': '1',
                'data-slider-step': '1',
                'data-slider-value': '1'}),
        }


class UserStreamSettingsForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserStreamSettingsForm, self).__init__(*args, **kwargs)
        if 'instance' in kwargs:
            self.fields['stream_vector_weight'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].stream_vector_weight)
            self.fields['stream_author_weight'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].stream_author_weight)
            self.fields['stream_journal_weight'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].stream_journal_weight)

    class Meta:
        model = UserSettings
        fields = (
                  'stream_vector_weight',
                  'stream_author_weight',
                  'stream_journal_weight',
                  )
        widgets = {
            'stream_vector_weight': forms.TextInput(attrs={
                'data-slider-id': 'stream_vector_weight_slider',
                'type': 'text',
                'data-slider-min': '0',
                'data-slider-max': '1',
                'data-slider-step': '.01',
                'data-slider-value': '1'}),
            'stream_author_weight': forms.TextInput(attrs={
                'data-slider-id': 'stream_author_weight_slider',
                'type': 'text',
                'data-slider-min': '0',
                'data-slider-max': '1',
                'data-slider-step': '.01',
                'data-slider-value': '1'}),
            'stream_journal_weight': forms.TextInput(attrs={
                'data-slider-id': 'stream_journal_weight_slider',
                'type': 'text',
                'data-slider-min': '0',
                'data-slider-max': '1',
                'data-slider-step': '.01',
                'data-slider-value': '1'}),
        }


class UserTrendSettingsForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserTrendSettingsForm, self).__init__(*args, **kwargs)
        if 'instance' in kwargs:
            self.fields['trend_doc_weight'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].trend_doc_weight)
            self.fields['trend_altmetric_weight'].widget.attrs['data-slider-value'] = \
                '{0:.2f}'.format(kwargs['instance'].trend_altmetric_weight)

    class Meta:
        model = UserSettings
        fields = (
            'trend_doc_weight',
            'trend_altmetric_weight',
        )
        widgets = {
            'trend_doc_weight': forms.TextInput(attrs={
                'data-slider-id': 'trend_doc_weight_slider',
                'type': 'text',
                'data-slider-min': '0',
                'data-slider-max': '1',
                'data-slider-step': '.01',
                'data-slider-value': '1'}),
            'trend_altmetric_weight': forms.TextInput(attrs={
                'data-slider-id': 'trend_altmetric_weight_slider',
                'type': 'text',
                'data-slider-min': '0',
                'data-slider-max': '1',
                'data-slider-step': '.01',
                'data-slider-value': '1'}),
        }


class UserEmailDigestSettingsForm(forms.ModelForm):

    class Meta:
        model = UserSettings
        fields = ('email_digest_frequency',
                  )
        widgets = {
            'email_digest_frequency': forms.Select(attrs={'class': 'form-control'}),
        }
