# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.contrib.auth import get_user_model

from .models import Thread

User = get_user_model()


class ThreadCreateForm(forms.ModelForm):

    class Meta:
        model = Thread
        fields = ('type', 'title', 'paper', 'author')
        exclude = ('author',)

    def __init__(self, *args, **kwargs):
        super(ThreadCreateForm, self).__init__(*args, **kwargs)
        self.fields['title'].required = True


class ThreadUpdateForm(ThreadCreateForm):

    class Meta(ThreadCreateForm.Meta):
        fields = ('title', 'content')

