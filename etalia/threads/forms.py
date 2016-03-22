# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import forms
from django.contrib.auth import get_user_model

from .models import Thread

User = get_user_model()


class ThreadForm(forms.ModelForm):

    class Meta:
        model = Thread
        field = ('type', 'title', 'paper', 'author')