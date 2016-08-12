# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import render, redirect
from django.views.decorators.cache import cache_page
from .tasks import failing_task


def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:my_feeds')
    else:
        context = {}
        return render(request, 'landing.html', context=context)


@cache_page(60 * 60)
def about(request):
    context = {}
    return render(request, 'pages/about.html', context=context)


@cache_page(60 * 60)
def terms_use(request):
    context = {}
    return render(request, 'pages/terms-use.html', context=context)


@cache_page(60 * 60)
def terms_privacy(request):
    context = {}
    return render(request, 'pages/terms-privacy.html', context=context)


@cache_page(60 * 60)
def support(request):
    context = {}
    return render(request, 'pages/support.html', context=context)


@cache_page(60 * 60)
def contact(request):
    context = {}
    return render(request, 'pages/contact.html', context=context)


@cache_page(60 * 60)
def help(request):
    context = {}
    return render(request, 'pages/help.html', context=context)


def test(request):
    return render(request, 'test.html', {})


def test_failing_task(request):
    failing_task.delay()
    return redirect('feeds:my_feeds')




