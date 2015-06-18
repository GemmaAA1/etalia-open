from django.shortcuts import render, redirect
from django.core.urlresolvers import reverse

def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:home')
    else:
        return render(request, 'landing.html')

def test(request):
    return render(request, 'test.html', {})