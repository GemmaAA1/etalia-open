# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import render, redirect
from django.db.models import Q
from django.core.urlresolvers import reverse

def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:main')
    else:
        return render(request, 'landing.html')

def test(request):
    return render(request, 'test.html', {})