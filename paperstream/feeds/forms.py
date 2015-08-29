# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.core.exceptions import ValidationError
from .models import UserFeed

class CreateUserFeedForm(forms.ModelForm):

    class Meta:
        model = UserFeed
        fields = ('name', )
        widgets = {
            'name':
                forms.TextInput(attrs={'class': 'form-control input-md ',
                                                'placeholder': 'Name'}),
        }

    def __init__(self, *args, **kwargs):
        super(CreateUserFeedForm, self).__init__(*args, **kwargs)
        self.fields['name'].initial = ''
