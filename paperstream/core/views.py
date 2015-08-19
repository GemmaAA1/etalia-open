from django.shortcuts import render, redirect
from django.db.models import Q
from django.core.urlresolvers import reverse

def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:feed')
    else:
        return render(request, 'landing.html')

def test(request):
    return render(request, 'test.html', {})