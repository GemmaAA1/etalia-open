# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth import get_user_model

from .models import Thread
from .constant import THREAD_QUESTION, THREAD_PAPER

User = get_user_model()


class ThreadCreateForm(forms.ModelForm):

    class Meta:
        model = Thread
        fields = ('type', 'title', 'paper')

    def __init__(self, *args, **kwargs):
        super(ThreadCreateForm, self).__init__(*args, **kwargs)
        self.fields['paper'].queryset = self.initial['papers']

    def clean_paper(self):
        cleaned_data = self.cleaned_data
        paper = cleaned_data.get('paper')
        type = cleaned_data.get('type')
        if type == THREAD_PAPER:    # paper is required
            if not paper:
                raise ValidationError('A publication is required for this thread type')
        return paper

    def clean(self):
        cleaned_data = self.cleaned_data
        # Check for unique_together
        if Thread.objects.filter(user_id=self.initial['user_id'],
                                 **cleaned_data).exists():
            if cleaned_data.get('type') == THREAD_QUESTION:
                raise ValidationError('A Thread with this title already exists for you')
            elif cleaned_data.get('type') == THREAD_PAPER:
                raise ValidationError('A Thread with this title and paper already exists for you')
            else:
                raise ValidationError('Unknown type of thread')
        else:
            return cleaned_data


class ThreadUpdateForm(ThreadCreateForm):

    class Meta(ThreadCreateForm.Meta):
        fields = ('title', 'content')

