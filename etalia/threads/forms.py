# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from .models import Thread, PubPeer, PubPeerComment


class ThreadForm(forms.ModelForm):

    class Meta:
        model = Thread
        fields = (
            'type',
            'user',
            'privacy',
            'paper',
            'title',
            'content'
        )

    def clean_title(self, data):
        title = self.cleaned_data['title']
        # if title is longer than max_length, truncate
        title_max = self.Meta.model._meta.get_field('title').max_length
        if len(title) > title_max:
            title = title[:title_max]
        return title


class PubPeerForm(forms.ModelForm):

    class Meta:
        model = PubPeer
        fields = (
            'thread',
            'doi',
            'link',
            'pubpeer_id',
        )


class PubPeerCommentForm(forms.ModelForm):

    class Meta:
        model = PubPeerComment
        fields = (
            'pubpeer',
            'body',
            'date',
            'pubpeercomment_id',
            'permalink',
            'rating',
            'user',
        )

