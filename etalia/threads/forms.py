# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.auth import get_user_model

from .models import Thread, ThreadPost, ThreadPostComment, ThreadMember
from .constant import THREAD_QUESTION, THREAD_PAPER

User = get_user_model()


class ThreadCreateForm(forms.ModelForm):
    class Meta:
        model = Thread
        fields = ('type', 'privacy', 'title', 'paper')

    def __init__(self, *args, **kwargs):
        # pop extra kwargs
        self.paper_qs = kwargs.pop('paper_qs')
        self.user_id = kwargs.pop('user_id')
        # init
        super(ThreadCreateForm, self).__init__(*args, **kwargs)
        # prepopulate
        self.fields['paper'].queryset = self.paper_qs

    def clean_paper(self):
        cleaned_data = self.cleaned_data
        paper = cleaned_data.get('paper')
        type = cleaned_data.get('type')
        if type == THREAD_PAPER:  # paper is required
            if not paper:
                raise ValidationError(
                    'A publication is required for this thread type')
        return paper

    def clean(self):
        cleaned_data = self.cleaned_data
        # Check coherence between 'type' and 'paper'
        if cleaned_data.get('type') == THREAD_QUESTION:
            if cleaned_data.get('paper'):
                raise ValidationError(
                    'Paper must be empty with {type} type thread. (and not displayed)'.format(
                        type=cleaned_data.get('type')))
        elif cleaned_data.get('type') == THREAD_PAPER:
            if not cleaned_data.get('paper'):
                raise ValidationError(
                    'Paper cannot be left empty with {type} type thread'.format(
                        type=cleaned_data.get('type')))

        # Check for unique_together
        if Thread.objects.filter(owner_id=self.user_id,
                                 **cleaned_data).exists():
            if cleaned_data.get('type') == THREAD_QUESTION:
                raise ValidationError(
                    'A Thread with this title already exists for you')
            elif cleaned_data.get('type') == THREAD_PAPER:
                raise ValidationError(
                    'A Thread with this title and paper already exists for you')
            else:
                raise ValidationError('Unknown type of thread')
        else:
            return cleaned_data


class ThreadUpdateForm(forms.ModelForm):
    class Meta(ThreadCreateForm.Meta):
        fields = ('title', 'content')

    def __init__(self, *args, **kwargs):
        self.instance = kwargs['instance']
        self.user_id = kwargs.pop('user_id')
        # init
        super(ThreadUpdateForm, self).__init__(*args, **kwargs)

    def clean(self):
        cleaned_data = self.cleaned_data
        # Check for unique_together
        if Thread.objects.filter(owner_id=self.instance.owner.id, **cleaned_data)\
                .exclude(pk=self.instance.id)\
                .exists():
            if self.instance.type == THREAD_QUESTION:
                raise ValidationError(
                    'A Thread with this title already exists for you')
            elif self.instance.type == THREAD_PAPER:
                raise ValidationError(
                    'A Thread with this title and paper already exists for you')
            else:
                raise ValidationError('Unknown type of thread')
        else:
            return cleaned_data


class ThreadPostForm(forms.ModelForm):
    class Meta:
        model = ThreadPost
        fields = ('content',)


class ThreadPostCommentForm(forms.ModelForm):
    class Meta:
        model = ThreadPostComment
        fields = ('content',)


class ThreadMemberForm(forms.ModelForm):

    class Meta:
        model = ThreadMember
        fields = ('thread', 'member')

