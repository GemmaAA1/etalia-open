# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import render, redirect
from django.views.decorators.cache import cache_page
from django.views.generic import DetailView, RedirectView
from django.utils.text import slugify
from etalia.press.models import Press
from .tasks import failing_task


def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:my_feeds')
    else:
        context = {}
        return render(request, 'landing.html', context=context)


def home_test(request):
    if request.user.is_authenticated():
        return redirect('feeds:my_feeds')
    else:
        context = {}
        return render(request, 'landing_test.html', context=context)


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
def press(request):
    context = {'press': Press.objects.all()}
    return render(request, 'pages/press.html', context=context)


class PressDetail(DetailView):

    template_name = 'pages/press_detail.html'
    model = Press

press_slug = PressDetail.as_view()


class PressDetailPk(RedirectView):
    """Redirect to slug press url"""

    permanent = True
    query_string = True
    pattern_name = 'core:press-slug'

    def get_redirect_url(self, *args, **kwargs):
        press = Press.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(press.title)
        return super(PressDetailPk, self).get_redirect_url(*args, **kwargs)

press_pk = PressDetailPk.as_view()


@cache_page(60 * 60)
def help(request):
    context = {}
    return render(request, 'pages/help.html', context=context)


def test(request):
    return render(request, 'test.html', {})


def test_failing_task(request):
    failing_task.delay()
    return redirect('feeds:my_feeds')




