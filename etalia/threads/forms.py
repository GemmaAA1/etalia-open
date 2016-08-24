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


class PubPeerForm(forms.ModelForm):

    class Meta:
        model = PubPeer
        fields = (
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

