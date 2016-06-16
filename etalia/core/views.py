# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import render, redirect
from django.db.models import Q
from django.utils import timezone
from django.conf import settings

from etalia.altmetric.models import AltmetricModel
from .tasks import failing_task


def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:my_feeds')
    else:
        # Get some trending altmetric matches
        d = timezone.datetime.now().date() - \
            timezone.timedelta(days=settings.LANDING_ACTIVE_PAPERS_TIME_IN_DAYS)
        ten_recent_most_active_paper = AltmetricModel.objects\
            .filter(Q(paper__date_ep__gt=d) |
                    (Q(paper__date_ep=None) & Q(paper__date_pp__gt=d)))\
            [:settings.LANDING_ACTIVE_PAPERS_NUMBER]
        papers = [p.paper for p in ten_recent_most_active_paper]
        context = {'active_papers': papers}
        return render(request, 'landing.html', context=context)


def about(request):
    context = {}
    return render(request, 'pages/about.html', context=context)


def terms(request):
    context = {}
    return render(request, 'pages/terms.html', context=context)


def support(request):
    context = {}
    return render(request, 'pages/support.html', context=context)


def help(request):
    context = {}
    return render(request, 'pages/help.html', context=context)


def test(request):
    return render(request, 'test.html', {})


def test_failing_task(request):
    failing_task.delay()
    return redirect('feeds:stream')




